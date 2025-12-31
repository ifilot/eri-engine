#include <eri/one_electron/overlap.h>
#include <eri/one_electron/kinetic.h>
#include <eri/one_electron/nuclear.h>

namespace eri::ops {

// defaults to Huzinaga
struct Overlap {
    static double eval(const basis::CGF& a, const basis::CGF& b) {
        return eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Huzinaga);
    }
};

struct OverlapHuzinaga {
    static double eval(const basis::CGF& a, const basis::CGF& b) {
        return eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Huzinaga);
    }
};

struct OverlapHellsing {
    static double eval(const basis::CGF& a, const basis::CGF& b) {
        return eri::one_electron::overlap(a, b, eri::enums::OverlapMethod::Hellsing);
    }
};

// defaults to Huzinaga
struct Kinetic {
    static double eval(const basis::CGF& a, const basis::CGF& b) {
        return eri::one_electron::kinetic(a, b, eri::enums::KineticMethod::Huzinaga);
    }
};

struct KineticHuzinaga {
    static double eval(const basis::CGF& a, const basis::CGF& b) {
        return eri::one_electron::kinetic(a, b, eri::enums::KineticMethod::Huzinaga);
    }
};

struct KineticHellsing {
    static double eval(const basis::CGF& a, const basis::CGF& b) {
        return eri::one_electron::kinetic(a, b, eri::enums::KineticMethod::Hellsing);
    }
};

// defaults to Huzinaga
struct Nuclear {
    static double eval(const basis::CGF& a, const basis::CGF& b, const std::array<double, 3>& C) {
        return eri::one_electron::nuclear(a, b, C, eri::enums::NuclearMethod::Huzinaga);
    }
};

struct NuclearHuzinaga {
    static double eval(const basis::CGF& a, const basis::CGF& b, const std::array<double, 3>& C) {
        return eri::one_electron::nuclear(a, b, C, eri::enums::NuclearMethod::Huzinaga);
    }
};

struct NuclearHellsing {
    static double eval(const basis::CGF& a, const basis::CGF& b, const std::array<double, 3>& C) {
        return eri::one_electron::nuclear(a, b, C, eri::enums::NuclearMethod::Hellsing);
    }
};

} // namespace eri::ops
