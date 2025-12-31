#include <eri/one_electron/overlap.h>

namespace eri::one_electron {

// forward declarations
double overlap_huzinaga(const eri::basis::CGF& a, const eri::basis::CGF& b, const eri::basis::CGF& c, const eri::basis::CGF& d);
double overlap_hellsing(const eri::basis::CGF& a, const eri::basis::CGF& b, const eri::basis::CGF& c, const eri::basis::CGF& d);

double overlap(const eri::basis::CGF& a, const eri::basis::CGF& b, , const eri::basis::CGF& c, const eri::basis::CGF& d, eri::enums::ERIMethod method) {
    // Dispatch lives here for now (we'll forward to Hellsing once implemented)
    switch (method) {
        case eri::enums::ERIMethod::Huzinaga:
            return eri_huzinaga(a, b);
        case eri::enums::ERIMethod::Hellsing:
            return eri_hellsing(a, b);
        default:
            __builtin_unreachable();
    }
}

} // namespace: eri::one_electron