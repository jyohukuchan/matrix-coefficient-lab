"""Small exact, executable counterparts of selected all-fields proof steps.

This package does not construct the schemes promised by the 9/4 theorem.
"""

from .fields import FiniteField
from .polynomials import Polynomial
from .schemes import BilinearScheme, descend, naive_scheme, strassen_scheme
from .degenerations import PolynomialDegeneration

__all__ = ["FiniteField", "Polynomial", "PolynomialDegeneration", "BilinearScheme",
           "descend", "naive_scheme", "strassen_scheme"]
