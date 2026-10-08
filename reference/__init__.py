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

__all__ = ["FiniteField", "Polynomial", "PolynomialDegeneration", "BilinearScheme",
           "FiniteTensor", "TensorScheme", "TensorDegeneration", "matrix_multiplication_tensor",
           "ThreeSectorConstruction", "convolution_tensor", "convolution_scheme", "lagrange_basis",
           "LinearMap", "shared_first_tensor", "tag_shared_scheme", "FiniteSeparation",
           "descend", "naive_scheme", "strassen_scheme"]
