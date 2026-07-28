#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

struct SimplePeakWindowResult {
	int peak_pos = -1;
	int left = -1;
	int right = -1;
	float height = 0.0f;
	bool ok = false;
};

struct SimpleTimeResult {
	bool ok = false;
	float sample_pos = 0.0f;
	float time_ns = 0.0f;
	bool used_fallback = false;
};

struct SimpleLogNormalFitResult {
	bool ok = false;
	float sample_pos = 0.0f;
	float time_ns = 0.0f;
	float chi2_per_dof = std::numeric_limits<float>::infinity();
	float A = 0.0f;
	float t0 = 0.0f;
	float mu = 0.0f;
	float sigma = 0.0f;
	std::string fail_reason;
	std::string equation;
};

inline float sampleToNs(float sample_pos, float clk_mhz)
{
	if (clk_mhz <= 0.0f) return sample_pos * 4.0f;
	return sample_pos * 1000.0f / clk_mhz;
}

inline float logNormalPulseValue(float t, float ped_mean,
							 float A, float t0, float mu, float sigma)
{
	if (!(sigma > 0.0f) || t <= t0) return ped_mean;
	const float dt = t - t0;
	const float z = (std::log(dt) - mu) / sigma;
	return ped_mean + A * std::exp(-0.5f * z * z);
}

inline std::string buildLogNormalEquation(float ped_mean,
							  const SimpleLogNormalFitResult &fit)
{
	if (!fit.ok) return "";
	return "y(t)= " + std::to_string(ped_mean)
		+ " + " + std::to_string(fit.A)
		+ " * exp(-0.5*((ln(t-" + std::to_string(fit.t0)
		+ ")-" + std::to_string(fit.mu)
		+ ")/" + std::to_string(fit.sigma) + ")^2), t>t0";
}

inline float evalLogNormalFitAtSample(float sample_idx, float ped_mean,
							   const SimpleLogNormalFitResult &fit)
{
	if (!fit.ok) return ped_mean;
	return logNormalPulseValue(sample_idx, ped_mean, fit.A, fit.t0, fit.mu, fit.sigma);
}

// Single-clean-peak helper:
// 1) find the maximum sample as peak (raw sample position),
// 2) compute ped-sub peak height and tail_cut = break_ratio * height,
// 3) walk left/right from peak using the same below-threshold logic as findPeaks.
inline SimplePeakWindowResult findSinglePeakWindowByRatio(const float waveform[100],
														  float ped_mean,
														  float ped_rms,
									  float break_ratio,
									  int smooth_order = 1)
{
	SimplePeakWindowResult out;
	if (waveform == nullptr) return out;

	if (break_ratio < 0.0f) break_ratio = 0.0f;
	if (break_ratio > 1.0f) break_ratio = 1.0f;

	constexpr int N = 100;
	int peak_pos = 0;
	float peak_val = waveform[0];
	for (int i = 1; i < N; ++i) {
		if (waveform[i] > peak_val) {
			peak_val = waveform[i];
			peak_pos = i;
		}
	}

	// Small correction from the original findPeaks:
	// search around the (smoothed) peak and pick the highest raw sample.
	int flat_left = peak_pos;
	while (flat_left > 0 && waveform[flat_left - 1] == waveform[peak_pos])
		--flat_left;
	int flat_right = peak_pos;
	while (flat_right < N - 1 && waveform[flat_right + 1] == waveform[peak_pos])
		++flat_right;

	const int search = std::max(1, smooth_order) + (flat_right - flat_left) / 2;
	int raw_pos = peak_pos;
	float raw_height = waveform[peak_pos] - ped_mean;
	for (int j = 1; j <= search; ++j) {
		if (peak_pos - j >= 0) {
			const float h = waveform[peak_pos - j] - ped_mean;
			if (h > raw_height) { raw_height = h; raw_pos = peak_pos - j; }
		}
		if (peak_pos + j < N) {
			const float h = waveform[peak_pos + j] - ped_mean;
			if (h > raw_height) { raw_height = h; raw_pos = peak_pos + j; }
		}
	}

	const float peak_height = raw_height;
	out.peak_pos = raw_pos;
	out.height = peak_height;
	if (peak_height <= 0.0f) return out;

	const float tail_cut = break_ratio * peak_height;
	auto is_below = [&](float v) {
		return v < tail_cut || v < ped_rms || v * peak_height < 0.0f;
	};

	int left = raw_pos;
	while (left > 0) {
		const float v = waveform[left - 1] - ped_mean;
		if (!is_below(v)) {
			--left;
		} else {
			break;
		}
	}

	int right = raw_pos;
	while (right < N - 1) {
		const float v = waveform[right + 1] - ped_mean;
		if (!is_below(v)) {
			++right;
		} else {
			break;
		}
	}

	out.left = left;
	out.right = right;
	out.ok = (out.left <= raw_pos && raw_pos <= out.right);
	return out;
}

// Time method 1: 3-point quadratic interpolation around peak sample.
inline SimpleTimeResult findTimeQuadraticNs(const float samples[100],
								int peak_pos,
								float clk_mhz)
{
	SimpleTimeResult out;
	if (samples == nullptr) return out;
	if (peak_pos < 0 || peak_pos >= 100) return out;

	float t_subsample = 0.0f;
	if (peak_pos > 0 && peak_pos < 99) {
		const float h_minus = samples[peak_pos - 1];
		const float h_zero  = samples[peak_pos];
		const float h_plus  = samples[peak_pos + 1];
		const float denom = h_minus - 2.0f * h_zero + h_plus;
		if (denom < -1.0e-3f) {
			const float delta = 0.5f * (h_minus - h_plus) / denom;
			t_subsample = std::max(-1.0f, std::min(1.0f, delta));
		}
	}

	out.sample_pos = static_cast<float>(peak_pos) + t_subsample;
	out.time_ns = sampleToNs(out.sample_pos, clk_mhz);
	out.ok = true;
	out.used_fallback = false;
	return out;
}

// Time method 2: linear CFD crossing on rising edge; fallback to quadratic.
inline SimpleTimeResult findTimeLinearCFDNs(const float samples[100],
								float ped_mean,
								float ped_rms,
								int left_bound,
								int peak_pos,
								float clk_mhz,
								float cfd_fraction = 0.5f)
{
	SimpleTimeResult quad = findTimeQuadraticNs(samples, peak_pos, clk_mhz);
	quad.used_fallback = true;
	if (samples == nullptr) return quad;
	if (peak_pos <= 0 || peak_pos >= 100) return quad;
	if (!(cfd_fraction > 0.0f) || !(cfd_fraction < 1.0f)) return quad;

	const int left = std::max(0, std::min(left_bound, peak_pos - 1));
	const float raw_height = samples[peak_pos] - ped_mean;
	if (!(raw_height > 0.0f)) return quad;

	const float v_thr = cfd_fraction * raw_height;
	bool cfd_ok = false;
	float pickoff = quad.sample_pos;

	for (int j = left-1; j < peak_pos; ++j) {
		const float v0 = samples[j] - ped_mean;
		const float v1 = samples[j + 1] - ped_mean;
		if (v0 < v_thr && v1 >= v_thr) {
			const float dv = v1 - v0;
			if (dv > 6.0f * ped_rms) {
				const float frac = (v_thr - v0) / dv;
				pickoff = j + std::max(0.0f, std::min(1.0f, frac));
				cfd_ok = true;
			}
		}
	}

	SimpleTimeResult out;
	out.ok = cfd_ok || quad.ok;
	out.sample_pos = cfd_ok ? pickoff : quad.sample_pos;
	out.time_ns = sampleToNs(out.sample_pos, clk_mhz);
	out.used_fallback = !cfd_ok;
	return out;
}

inline float chol4(float M[16], float L[16])
{
	float min_pivot_sq = std::numeric_limits<float>::infinity();
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j <= i; ++j) {
			float sum = M[i * 4 + j];
			for (int k = 0; k < j; ++k)
				sum -= L[i * 4 + k] * L[j * 4 + k];
			if (i == j) {
				if (sum <= 0.0f) return -1.0f;
				if (sum < min_pivot_sq) min_pivot_sq = sum;
				L[i * 4 + j] = std::sqrt(sum);
			} else {
				L[i * 4 + j] = sum / L[j * 4 + j];
			}
		}
		for (int j = i + 1; j < 4; ++j) L[i * 4 + j] = 0.0f;
	}
	return min_pivot_sq;
}

inline void chol4Solve(const float L[16], const float b[4], float x[4])
{
	float y[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	for (int i = 0; i < 4; ++i) {
		float s = b[i];
		for (int k = 0; k < i; ++k) s -= L[i * 4 + k] * y[k];
		y[i] = s / L[i * 4 + i];
	}
	for (int i = 3; i >= 0; --i) {
		float s = y[i];
		for (int k = i + 1; k < 4; ++k) s -= L[k * 4 + i] * x[k];
		x[i] = s / L[i * 4 + i];
	}
}

// Time method 3: local log-normal CFD fit, plus fit equation output.
inline SimpleLogNormalFitResult findTimeLogNormalCFDNs(const float samples[100],
									float ped_mean,
									float ped_rms,
									int left_bound,
									int peak_pos,
									float clk_mhz,
									float cfd_fraction = 0.5f)
{
	SimpleLogNormalFitResult out;
	auto fail = [&](const char *reason) -> SimpleLogNormalFitResult {
		out.fail_reason = reason;
		return out;
	};
	if (samples == nullptr) return fail("null samples");
	if (peak_pos < 0 || peak_pos >= 100) return fail("peak_pos out of range");
	if (!(cfd_fraction > 0.0f) || !(cfd_fraction < 1.0f)) return fail("invalid cfd_fraction");

	const int fit_start = std::max(0, left_bound - 5);
	const int fit_stop  = std::min(99, peak_pos + 4);
	const int nfit = fit_stop - fit_start + 1;
	if (nfit < 8) return fail("fit window too short (nfit < 8)");

	const float raw_height = samples[peak_pos] - ped_mean;
	if (!(raw_height > 0.0f)) return fail("non-positive raw height");

	const float v_thr = cfd_fraction * raw_height;
	int cfd_left_sample = -1;
	for (int j = std::max(0, left_bound-2); j < peak_pos; ++j) {
		const float v0 = samples[j] - ped_mean;
		const float v1 = samples[j + 1] - ped_mean;
		if (v0 < v_thr && v1 >= v_thr) cfd_left_sample = j;
	}
	if (cfd_left_sample < 0) return fail("no CFD crossing found before peak");

	float t[100];
	float y[100];
	for (int i = 0; i < nfit; ++i) {
		const int idx = fit_start + i;
		t[i] = static_cast<float>(idx);
		y[i] = samples[idx];
	}

	const float peak_time = static_cast<float>(peak_pos);
	const float A_guess = std::max(samples[peak_pos] - ped_mean, 1.0f);

	float p[4] = {
		A_guess,
		std::max(0.0f, peak_time - 15.0f / 4.0f),
		1.5f,
		0.6f,
	};
	float p_lo[4] = {
		0.0f,
		std::max(0.0f, t[0] - 12.0f / 4.0f),
		0.3f,
		0.1f,
	};
	float p_hi[4] = {
		3.0f * A_guess,
		peak_time - 1.0e-3f,
		4.0f,
		1.0f,
	};
	for (int j = 0; j < 4; ++j)
		p[j] = std::clamp(p[j], p_lo[j], p_hi[j]);

	const float sigma_fit = std::max(ped_rms, 1.0f);
	const float inv_sigma2 = 1.0f / (sigma_fit * sigma_fit);
	const int dof = std::max(1, nfit - 4);

	auto eval_chi2 = [&](const float *params) -> float {
		if (!(params[0] > 0.0f) || !(params[3] > 0.0f) || !(params[1] < peak_time))
			return std::numeric_limits<float>::infinity();
		float chi2 = 0.0f;
		for (int i = 0; i < nfit; ++i) {
			const float r = y[i] - logNormalPulseValue(
				t[i], ped_mean, params[0], params[1], params[2], params[3]);
			chi2 += r * r;
		}
		return (chi2 * inv_sigma2) / static_cast<float>(dof);
	};

	float chi2 = eval_chi2(p);
	if (!std::isfinite(chi2)) return fail("initial chi2 is not finite");

	float p_best[4];
	for (int j = 0; j < 4; ++j) p_best[j] = p[j];
	float chi2_best = chi2;

	float lambda = 1.0e-3f;
	bool any_accepted = false;

	float r[100];
	float J[4 * 100];
	float A[16];
	float Aug[16];
	float L[16];
	float g[4], delta[4], p_new[4];

	for (int iter = 0; iter < 60; ++iter) {
		for (int i = 0; i < nfit; ++i)
			r[i] = y[i] - logNormalPulseValue(t[i], ped_mean, p[0], p[1], p[2], p[3]);

		const float h_A     = std::max(1.0e-3f * std::max(std::abs(p[0]), 1.0f), 1.0e-3f);
		const float h_t0    = 1.0e-3f;
		const float h_mu    = std::max(1.0e-3f * std::max(std::abs(p[2]), 1.0f), 1.0e-4f);
		const float h_sigma = std::max(1.0e-3f * std::max(std::abs(p[3]), 1.0f), 1.0e-4f);
		const float h[4] = {h_A, h_t0, h_mu, h_sigma};

		for (int j = 0; j < 4; ++j) {
			float p_fd[4];
			for (int k = 0; k < 4; ++k) p_fd[k] = p[k];
			p_fd[j] = std::min(p_hi[j], p[j] + h[j]);
			const float step = p_fd[j] - p[j];
			if (!(step > 0.0f)) {
				for (int i = 0; i < nfit; ++i) J[j * 100 + i] = 0.0f;
				continue;
			}
			for (int i = 0; i < nfit; ++i) {
				const float f0 = logNormalPulseValue(t[i], ped_mean, p[0],    p[1],    p[2],    p[3]);
				const float fp = logNormalPulseValue(t[i], ped_mean, p_fd[0], p_fd[1], p_fd[2], p_fd[3]);
				J[j * 100 + i] = -(fp - f0) / step;
			}
		}

		for (int j = 0; j < 4; ++j) {
			g[j] = 0.0f;
			for (int k = 0; k < 4; ++k) A[j * 4 + k] = 0.0f;
		}
		for (int i = 0; i < nfit; ++i) {
			for (int j1 = 0; j1 < 4; ++j1) {
				const float jj1 = J[j1 * 100 + i];
				g[j1] += jj1 * r[i] * inv_sigma2;
				for (int j2 = 0; j2 <= j1; ++j2) {
					const float jj2 = J[j2 * 100 + i];
					A[j1 * 4 + j2] += jj1 * jj2 * inv_sigma2;
				}
			}
		}
		for (int j1 = 0; j1 < 4; ++j1)
			for (int j2 = j1 + 1; j2 < 4; ++j2)
				A[j1 * 4 + j2] = A[j2 * 4 + j1];

		for (int j = 0; j < 16; ++j) Aug[j] = A[j];
		for (int j = 0; j < 4; ++j) Aug[j * 4 + j] += lambda;

		const float min_pivot_sq = chol4(Aug, L);
		if (min_pivot_sq < 0.0f) {
			lambda *= 10.0f;
			if (lambda > 1.0e10f) break;
			continue;
		}
		chol4Solve(L, g, delta);

		float dpar = 0.0f;
		for (int j = 0; j < 4; ++j) {
			p_new[j] = std::clamp(p[j] - delta[j], p_lo[j], p_hi[j]);
			const float scale = std::max(std::abs(p[j]), 1.0e-3f);
			dpar += std::abs(p_new[j] - p[j]) / scale;
		}

		const float chi2_new = eval_chi2(p_new);
		if (std::isfinite(chi2_new) && chi2_new < chi2_best) {
			chi2_best = chi2_new;
			for (int j = 0; j < 4; ++j) p_best[j] = p_new[j];
		}

		if (std::isfinite(chi2_new) && chi2_new < chi2) {
			for (int j = 0; j < 4; ++j) p[j] = p_new[j];
			chi2 = chi2_new;
			any_accepted = true;
			lambda = std::max(lambda / 10.0f, 1.0e-12f);
			if (dpar < 1.0e-4f) break;
		} else {
			lambda *= 10.0f;
			if (lambda > 1.0e10f) break;
		}
	}

	if (!any_accepted && !std::isfinite(chi2_best)) return fail("LM never accepted a step");

	const float root_term = std::sqrt(-2.0f * std::log(cfd_fraction));
	const float t_cfd_sample = p_best[1] + std::exp(p_best[2] - p_best[3] * root_term);
	const float t_peak_sample = p_best[1] + std::exp(p_best[2]);
	const float fit_lo_sample = static_cast<float>(fit_start);
	const float fit_hi_sample = static_cast<float>(peak_pos);
	const float cfd_lo_sample = static_cast<float>(cfd_left_sample);
	const float cfd_hi_sample = cfd_lo_sample + 1.0f;

	if (!std::isfinite(t_cfd_sample) || !std::isfinite(t_peak_sample)) return fail("final time is not finite");
	//if (!(t_cfd_sample >= fit_lo_sample && t_cfd_sample <= fit_hi_sample)) return fail("t_cfd outside fit window");
	//if (!(t_cfd_sample >= cfd_lo_sample && t_cfd_sample <= cfd_hi_sample)) return fail("t_cfd outside CFD bracket");
	//if (!(t_peak_sample >= fit_lo_sample && t_peak_sample <= static_cast<float>(fit_stop) + 1.0f)) return fail("t_peak outside fit window");

	out.ok = true;
	out.sample_pos = t_cfd_sample;
	out.time_ns = sampleToNs(t_cfd_sample, clk_mhz);
	out.chi2_per_dof = chi2_best;
	out.A = p_best[0];
	out.t0 = p_best[1];
	out.mu = p_best[2];
	out.sigma = p_best[3];
	out.equation = buildLogNormalEquation(ped_mean, out);
	return out;
}
