#include "../../include/eri/cgf.h"
#include "../../include/eri/one_electron/overlap.h"

namespace eri::one_electron {

// forward declarations
double overlap_huzinaga(const eri::basis::CGF& a, const eri::basis::CGF& b);
double overlap_hellsing(const eri::basis::CGF& a, const eri::basis::CGF& b);

double overlap(const eri::basis::CGF& a, const eri::basis::CGF& b, OverlapMethod method) {
    // Dispatch lives here for now (we'll forward to Hellsing once implemented)
    switch (method) {
        case OverlapMethod::Huzinaga:
            return overlap_huzinaga(a, b);
        case OverlapMethod::Hellsing:
            // placeholder until we add Hellsing implementation file below
            return overlap_huzinaga(a, b);
        default:
            __builtin_unreachable();
    }
}

} // namespace: eri::one_electron