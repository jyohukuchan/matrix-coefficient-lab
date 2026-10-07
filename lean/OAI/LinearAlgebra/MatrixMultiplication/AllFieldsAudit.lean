import OAI.LinearAlgebra.MatrixMultiplication.Main
import Mathlib.Data.ZMod.Basic

/-!
# Checks of the all-fields theorem and its specification

The generic example below has only `[Field F]`, in an arbitrary universe.
The small-characteristic examples also exercise the finite fields for which
the original fixed Fourier period or integer interpolation nodes can fail.
The cost example retains the upstream program correctness and cost model.
-/

namespace OAI.MatrixMultiplication.AllFieldsAudit

universe u

example (F : Type u) [Field F] : Arithmetic.omega F ≤ (9 : ℝ) / 4 :=
  omega_le_nine_quarters F

example : Arithmetic.omega (ZMod 2) ≤ (9 : ℝ) / 4 :=
  omega_le_nine_quarters (ZMod 2)

example : Arithmetic.omega (ZMod 3) ≤ (9 : ℝ) / 4 :=
  omega_le_nine_quarters (ZMod 3)

example : Arithmetic.omega (ZMod 5) ≤ (9 : ℝ) / 4 :=
  omega_le_nine_quarters (ZMod 5)

example : Arithmetic.omega ℚ ≤ (9 : ℝ) / 4 := omega_le_nine_quarters ℚ

example : Arithmetic.omega ℝ ≤ (9 : ℝ) / 4 := omega_le_nine_quarters ℝ

example : Arithmetic.omega ℂ ≤ (9 : ℝ) / 4 := complex_omega_le_nine_quarters

example (F : Type u) [Field F] (ε : ℝ) (hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∀ n : ℕ, 1 ≤ n →
      ∃ P : Arithmetic.MatrixAlgorithm F n n n,
        P.Correct ∧ (P.cost : ℝ) ≤ C * (n : ℝ) ^ ((9 : ℝ) / 4 + ε) :=
  AuxiliarySeparation.matrix_multiplication_cost_le F ε hε

#check @OAI.MatrixMultiplication.omega_le_nine_quarters
#print axioms OAI.MatrixMultiplication.omega_le_nine_quarters
#print axioms OAI.MatrixMultiplication.complex_omega_le_nine_quarters
#print axioms OAI.MatrixMultiplication.AuxiliarySeparation.matrix_multiplication_cost_le
#print axioms OAI.MatrixMultiplication.AuxiliarySeparation.exactRankExponent_le_algebraicClosure
#print axioms OAI.MatrixMultiplication.AuxiliarySeparation.exists_detecting_character
#print axioms OAI.MatrixMultiplication.Foundation.Tensor.RankAtMost.descend_algebraicClosure_powers

end OAI.MatrixMultiplication.AllFieldsAudit
