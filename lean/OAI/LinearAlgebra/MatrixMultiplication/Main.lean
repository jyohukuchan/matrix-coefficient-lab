/-
Upstream source and modification notice
Modified from OpenAI's matrix-multiplication sources to support arbitrary fields.
See UPSTREAM.md for provenance and changes; Apache-2.0 license retained in LICENSE.
-/

import OAI.LinearAlgebra.MatrixMultiplication.AuxiliarySeparation.Main
import OAI.LinearAlgebra.MatrixMultiplication.ComplexBounds.CW75ReorderedRectangular
import OAI.LinearAlgebra.MatrixMultiplication.Duality.ExponentBound
import OAI.LinearAlgebra.MatrixMultiplication.FieldConstruction.Conditional
import OAI.LinearAlgebra.MatrixMultiplication.Numerical.AllFieldCertificatesReplayIntegration

namespace OAI

/-! Arithmetic complexity bounds for square and rectangular matrix multiplication. -/

namespace MatrixMultiplication

/-- The square matrix multiplication exponent over every field is at most `9/4`. -/
theorem omega_le_nine_quarters (F : Type*) [Field F] :
    Arithmetic.omega F ≤ (9 : ℝ) / 4 :=
  AuxiliarySeparation.omega_le_nine_quarters F

/-- The complex square matrix multiplication exponent is at most `9/4`. -/
theorem complex_omega_le_nine_quarters :
    Arithmetic.omega ℂ ≤ (9 : ℝ) / 4 :=
  omega_le_nine_quarters ℂ

/-- The complex dual matrix multiplication exponent is greater than `0.465`. -/
theorem complex_alpha_gt_93_div_200 :
    (93 : ℝ) / 200 < Arithmetic.complexAlpha :=
  DualExponentBound.alpha_gt

/-- Multiplication of `n × n^0.709` by `n^0.709 × n` complex matrices has
arithmetic exponent less than `2.092`. -/
theorem complex_rectangular_omega_lt_523_div_250 :
    Arithmetic.rectangularOmega ℂ ((709 : ℝ) / 1000) < (523 : ℝ) / 250 :=
  CW75ReorderedRectangular.rectangularOmega_lt_target

end MatrixMultiplication

end OAI
