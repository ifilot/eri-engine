#include <eri/basis/basisset.h>

#include <fstream>
#include <stdexcept>

#include <json/json.h>

namespace eri::basis {

/**
 * @brief  Autogenerates angular momentum values
 * @note   
 * @param  l: angular momentum sum (0->s, 1->p, 2->d, etc)
 * @retval Angular momenta values
 */
static std::vector<std::array<int,3>>
cartesian_angular_momenta(int l)
{
    std::vector<std::array<int,3>> am;
    for (int lx = l; lx >= 0; --lx)
        for (int ly = l - lx; ly >= 0; --ly) {
            int lz = l - lx - ly;
            am.push_back({lx, ly, lz});
        }
    return am;
}

/**
 * @brief  Construct basis set from .json file
 * @note   Uses Basis Set Exchange standard
 * @param  filename: basis set file
 * @param  mol: Molecule
 * @retval None
 */
void BasisSet::load_from_bse_json(
    const std::string& filename,
    const eri::chem::Molecule& mol
) {
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("Failed to open basis file:" + filename);

    Json::Value root;
    in >> root;

    const Json::Value& elements = root["elements"];

    for (const auto& atom : mol.atoms()) {
        const std::string Zstr = std::to_string(atom.Z);
        const Json::Value& shells = elements[Zstr]["electron_shells"];

        for (const auto& shell : shells) {

            const Json::Value& am_list = shell["angular_momentum"];
            const Json::Value& exps_json = shell["exponents"];
            const Json::Value& coefs_json = shell["coefficients"];

            // Parse exponents once per shell
            std::vector<double> exponents;
            for (const auto& e : exps_json)
                exponents.push_back(std::stod(e.asString()));

            for (Json::ArrayIndex i = 0; i < am_list.size(); ++i) {
                int l = am_list[i].asInt();

                // Parse contraction coefficients for this angular momentum
                std::vector<double> coefficients;
                for (const auto& c : coefs_json[i])
                    coefficients.push_back(std::stod(c.asString()));

                auto cart_am = cartesian_angular_momenta(l);

                for (const auto& [lx, ly, lz] : cart_am) {
                    cgfs_.emplace_back(
                        lx, ly, lz,
                        atom.position,
                        exponents,
                        coefficients
                    );
                }
            }
        }
    }
}

} // namespace eri::basis
