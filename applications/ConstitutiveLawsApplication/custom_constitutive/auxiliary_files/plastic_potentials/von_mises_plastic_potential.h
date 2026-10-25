// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//

#pragma once
#include <cmath>

// System includes

// Project includes
#include "generic_plastic_potential.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @class VonMisesPlasticPotential
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a plastic potential following the theory of Von Mises
 * @details If the plastic potential is of vonMises (cylinder) type, on can see that the plastic strain increment tensor is in principle the scaled deviatoric stress tensor, hence principal directions coincide. When the yield and plastic potential surfaces are plotted in principal stress space the resulting surface will be a circular cylinder for Von-Mises. This means that both yield and strength are dependent on intermediate principal stress, sigma_2
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <SizeType TVoigtSize = 6>
class VonMisesPlasticPotential
{
public:
    ///@name Type Definitions
    ///@{

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// The define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// Counted pointer of VonMisesPlasticPotential
    KRATOS_CLASS_POINTER_DEFINITION(VonMisesPlasticPotential);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Copy constructor
    VonMisesPlasticPotential(VonMisesPlasticPotential const &rOther)
    {
    }

    /// Assignment operator
    VonMisesPlasticPotential &operator=(VonMisesPlasticPotential const &rOther)
    {
        return *this;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief This  script  calculates  the derivatives  of the plastic potential
    according   to   NAYAK-ZIENKIEWICZ   paper International
    journal for numerical methods in engineering vol 113-135 1972.
    As:            DF/DS = c1*V1 + c2*V2 + c3*V3
     * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
     * @param rDeviator The deviatoric part of the stress vector
     * @param J2 The second invariant of the Deviator
     * @param rGFlux The derivative of the plastic potential
     * @param rValues Parameters of the constitutive law
     */
    static void CalculatePlasticPotentialDerivative(
        const array_1d<double, VoigtSize>& rPredictiveStressVector,
        const array_1d<double, VoigtSize>& rDeviator,
        const double J2,
        array_1d<double, VoigtSize>& rGFlux,
        ConstitutiveLaw::Parameters& rValues
        )
    {
        array_1d<double, VoigtSize> first_vector, second_vector, third_vector;
        AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateSecondVector(rDeviator, J2, second_vector);
        const double c2 = std::sqrt(3.0);

        noalias(rGFlux) = c2 * second_vector;
    }

    /** @brief Opt-in gradient and Hessian in stress/engineering-strain Voigt.
     * Requires q from this stress, q > 0. No absolute stress-unit cutoff.
     * Existing potential methods and callers are intentionally unchanged.
     */
    static bool CalculatePlasticPotentialDerivativeAndJacobian(
        const array_1d<double, VoigtSize>& rStress, double Q,
        array_1d<double, VoigtSize>& rG, Matrix& rDG)
    {
        if (VoigtSize != 6 || !std::isfinite(Q) || Q <= 0.0) return false;
        const double mean = rStress[0] / 3.0 + rStress[1] / 3.0 + rStress[2] / 3.0;
        for (SizeType i = 0; i < VoigtSize; ++i) {
            const double deviator = rStress[i] - (i < 3 ? mean : 0.0);
            rG[i] = (i < 3 ? 1.5 : 3.0) * (deviator / Q);
            if (!std::isfinite(rG[i])) return false;
        }
        rDG.resize(VoigtSize, VoigtSize, false);
        for (SizeType i = 0; i < VoigtSize; ++i) {
            for (SizeType j = 0; j < VoigtSize; ++j) {
                const double metric_projector = i < 3 && j < 3
                    ? (i == j ? 2.0 / 3.0 : -1.0 / 3.0)
                    : (i == j ? 2.0 : 0.0);
                rDG(i,j) = (1.5 * metric_projector - rG[i] * rG[j]) / Q;
                if (!std::isfinite(rDG(i,j))) return false;
            }
        }
        return true;
    }

    /**
     * @brief This method defines the check to be performed in the plastic potential
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        return 0;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class GenericYieldSurface

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
