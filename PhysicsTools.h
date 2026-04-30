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