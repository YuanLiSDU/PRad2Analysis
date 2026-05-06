// Physical constants
static constexpr float M_PROTON  = 938.272f;   // MeV
static constexpr float M_ELECTRON = 0.511f;    // MeV
static constexpr float DEG2RAD = 3.14159265f / 180.f;

// two simple data structures used in physics analysis
struct GEMHit {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    uint8_t det_id = 5; // 0-3 for GEM1-GEM4
};

struct HCHit {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    float energy = 0.f;
    uint16_t center_id = 0; // index of central block
    uint32_t flag = -1;

    HCHit() = default;
    HCHit(float x_, float y_, float z_, float e_) : x(x_), y(y_), z(z_), energy(e_) {}
};

//data structure for storing reconstructed Moller events used for analysis
struct DataPoint
{
    float x;
    float y;
    float z;
    float E;

    DataPoint() {};
    DataPoint(float xi, float yi, float zi, float Ei) : x(xi), y(yi), z(zi), E(Ei) {};
};
typedef std::pair<DataPoint, DataPoint> MollerEvent;
typedef std::vector<MollerEvent> MollerData;

std::array<float, 2> GetMollerCenter(const MollerEvent &event1, const MollerEvent &event2)
{
    float x1[2], y1[2];
    float x2[2], y2[2];

    x1[0] = event1.first.x;  y1[0] = event1.first.y;
    x1[1] = event1.second.x; y1[1] = event1.second.y;
    x2[0] = event2.first.x;  y2[0] = event2.first.y;
    x2[1] = event2.second.x; y2[1] = event2.second.y;

    //two lines: y = ax + b, y = cx + d
    float dx1 = x1[0] - x1[1];
    float dx2 = x2[0] - x2[1];
    if (std::abs(dx1) < 1e-6f || std::abs(dx2) < 1e-6f)
        return {0.f, 0.f};  // vertical line — degenerate

    float a = (y1[0] - y1[1]) / dx1;
    float b = y1[0] - a * x1[0];
    float c = (y2[0] - y2[1]) / dx2;
    float d = y2[0] - c * x2[0];

    if (std::abs(a - c) < 1e-6f)
        return {0.f, 0.f};  // parallel lines — no intersection

    float x_cross = (d - b) / (a - c);
    float y_cross = a * x_cross + b;

    return {x_cross, y_cross};

}

float GetMollerZdistance(const MollerEvent &event, float Ebeam)
{
    float R1 = sqrt(event.first.x*event.first.x + event.first.y*event.first.y);
    float R2 = sqrt(event.second.x*event.second.x + event.second.y*event.second.y);
    float z = sqrt( (Ebeam + M_ELECTRON) * R1 * R2 / (2.*M_ELECTRON) );
    return z;
}

float GetPhiAngle(float x, float y)
{
    // atan2 handles all quadrants and x==0 correctly
    float phi = std::atan2(y, x) * 180.f / static_cast<float>(TMath::Pi());
    if (phi < 0) phi += 360.f;
    return phi;
}

float GetMollerPhiDiff(const MollerEvent &event1)
{
    // Calculate the azimuthal angle difference (phi) for a Moller event
    float x1 = event1.first.x, y1 = event1.first.y;
    float x2 = event1.second.x, y2 = event1.second.y;
    float phi1 = GetPhiAngle(x1, y1);
    float phi2 = GetPhiAngle(x2, y2);
    float phi_diff = fabs(phi1 - phi2) - 180.f; // Expecting back-to-back, so difference should be around 180 degrees
    return phi_diff;
}

float ExpectedEnergy(float theta_deg, float EBeam, const std::string &type)
{
    float theta = theta_deg * DEG2RAD;
    float cos_t = std::cos(theta);
    float sin_t = std::sin(theta);

    if (type == "ep") {
        // elastic e-p: E' = E * M / (M + E*(1 - cos_t))
        // where M = proton mass
        float expectE = EBeam * M_PROTON / (M_PROTON + EBeam * (1.f - cos_t));
        return expectE;
    }
    if (type == "ee") {
        // Moller scattering: exact lab-frame formula from 4-momentum conservation
        // E' = m * [(gamma+1) + (gamma-1)*cos^2(theta)] / [(gamma+1) - (gamma-1)*cos^2(theta)]
        float gamma = EBeam / M_ELECTRON;
        float num = (gamma + 1.f) + (gamma - 1.f) * cos_t * cos_t;
        float den = (gamma + 1.f) - (gamma - 1.f) * cos_t * cos_t;
        if (den <= 0) return 0.f;
        float expectE = M_ELECTRON * num / den;
        return expectE;
    }
    return 0.f;
}

bool isMott(float energy, float EBeam, float sigma_percent)
{
    return fabs(energy - EBeam) < 3.f * sigma_percent * EBeam / sqrt(EBeam/1000.f);
}

bool isMoller_kinematic(float theta1, float energy1, float theta2, float energy2, float EBeam, float sigma_percent)
{
    float expectE1 = ExpectedEnergy(theta1, EBeam, "ee");
    float expectE2 = ExpectedEnergy(theta2, EBeam, "ee");

    bool E_sum = false, E1_ok = false, E2_ok = false, phi_ok = false;

    if(fabs(energy1 + energy2 - EBeam) < 4.f * sigma_percent * EBeam / sqrt(EBeam/1000.f)) 
        E_sum = true;
    if(fabs(energy1 - expectE1) < 5.f * expectE1 * sigma_percent / sqrt(expectE1/1000.f)) 
        E1_ok = true;
    if(fabs(energy2 - expectE2) < 5.f * expectE2 * sigma_percent / sqrt(expectE2/1000.f)) 
        E2_ok = true;

    return E_sum && E1_ok && E2_ok;
}