"""Small exact, executable counterparts of selected all-fields proof steps.

This package does not construct the schemes promised by the 9/4 theorem.
"""

from .fields import FiniteField
from .polynomials import Polynomial
from .schemes import BilinearScheme, descend, naive_scheme, strassen_scheme
from .degenerations import PolynomialDegeneration
from .tensors import FiniteTensor, matrix_multiplication_tensor, shared_first_tensor
from .tensor_schemes import TensorScheme
from .tensor_degenerations import TensorDegeneration
from .sectors import ThreeSectorConstruction, convolution_tensor
from .convolution import convolution_scheme, lagrange_basis
from .maps import LinearMap
from .tagging import tag_shared_scheme
from .separation import FiniteSeparation
from .extraction import (ProofMatrixConstruction, branch_maps, extract_branch,
                         factor_singleton_leg, proof_matrix_pipeline,
                         simplify_scheme, square_from_boundary)
from .catalysts import (CatalyticCertificate, ConditionalCostPlan,
                       plan_catalytic_parameters, plan_conditional_costs)
from .catalyst_compiler import (CompiledAbsorption, CompiledGainAbsorption, CompiledRestriction,
                               CompilerLimits, CompilerResourceLimit, PoweredCatalyst,
                               PoweredGainCatalyst, compile_absorption, compile_gain_absorption,
                               compile_gain_powered_catalyst, compile_powered_catalyst,
                               coordinate_catalyst, unit_from_nonzero)
from .witnessed_constraints import (DetectorConstraint, DualSearchResult, IntegerStateCertificate,
                                    StateAssemblyLimit, StateAssemblyLimits, WitnessedOrder,
                                    WitnessedStateSystem, assemble_catalyst, clear_rational_dual,
                                    find_nonnegative_dual, rank_upper_constraint)
from .constraint_generators import (WitnessedOrderExport, export_sector_order, lift_witnessed_order)
from .projective_convolution import projective_convolution_scheme
from .scalar_summands import eliminate_scalar_summands
from .gain_planning import (GainCompilationResult, GainCostPlan, compile_gain_target,
                            plan_gain_catalytic_parameters)
from .rectangular_restrictions import share_first_axis, share_matrix_left_operand
from .finite_type_constraints import (ExactTypeExport, export_exact_type_order,
                                     export_finite_separation_order)
from .determinant_filtration import DeterminantFiltration

__all__ = ["FiniteField", "Polynomial", "PolynomialDegeneration", "BilinearScheme",
           "FiniteTensor", "TensorScheme", "TensorDegeneration", "matrix_multiplication_tensor",
           "ThreeSectorConstruction", "convolution_tensor", "convolution_scheme", "lagrange_basis",
           "LinearMap", "shared_first_tensor", "tag_shared_scheme", "FiniteSeparation",
           "ProofMatrixConstruction", "branch_maps", "extract_branch", "factor_singleton_leg",
           "proof_matrix_pipeline", "simplify_scheme", "square_from_boundary",
           "CatalyticCertificate", "ConditionalCostPlan", "plan_catalytic_parameters",
           "plan_conditional_costs",
           "CompiledAbsorption", "CompiledGainAbsorption", "CompiledRestriction", "CompilerLimits",
           "CompilerResourceLimit", "PoweredCatalyst", "PoweredGainCatalyst", "compile_absorption",
           "compile_gain_absorption", "compile_gain_powered_catalyst", "compile_powered_catalyst", "coordinate_catalyst",
           "unit_from_nonzero", "DetectorConstraint", "DualSearchResult", "IntegerStateCertificate",
           "StateAssemblyLimit", "StateAssemblyLimits", "WitnessedOrder", "WitnessedStateSystem",
           "assemble_catalyst", "clear_rational_dual", "find_nonnegative_dual", "rank_upper_constraint",
           "WitnessedOrderExport", "export_sector_order", "lift_witnessed_order",
           "projective_convolution_scheme", "eliminate_scalar_summands", "share_first_axis",
           "GainCompilationResult", "GainCostPlan", "compile_gain_target", "plan_gain_catalytic_parameters",
           "ExactTypeExport", "export_exact_type_order", "export_finite_separation_order",
           "DeterminantFiltration",
           "share_matrix_left_operand",
           "descend", "naive_scheme", "strassen_scheme"]
