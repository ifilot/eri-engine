#pragma once

namespace eri::two_electron {

constexpr std::size_t pair_index(std::size_t i, std::size_t j) noexcept {
    return i * (i + 1) / 2 + j;
}

constexpr std::size_t eri_index(std::size_t i,
                                std::size_t j,
                                std::size_t k,
                                std::size_t l) noexcept
{
    // Order within pairs
    if (i < j) std::swap(i, j);
    if (k < l) std::swap(k, l);

    const std::size_t ij = pair_index(i, j);
    const std::size_t kl = pair_index(k, l);

    // Order bra/ket pairs
    const std::size_t p = ij >= kl ? ij : kl;
    const std::size_t q = ij >= kl ? kl : ij;

    return pair_index(p, q);
}

} // namespace eri::two_electron
