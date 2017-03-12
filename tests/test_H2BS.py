import unittest
import numpy as np
from spec_func import AssociatedLegendre


class AsLegTest(unittest.TestCase):
    def test_zero_order_asp(self):
        x = np.linspace(-3.0, 4.0)
        self.assertEqual(AssociatedLegendre(0, 0, x)[0], 1.0)
        self.assertEqual(AssociatedLegendre(0, 0, x)[x == 0], 1.0)
        self.assertEqual(AssociatedLegendre(0, 0, x)[-1], 1.0)

    def test_first_order_asp(self):
        x = np.linspace(-3.0, 4.0)
        self.assertEqual(AssociatedLegendre(1, 0, x)[x == -3.0], -3.0)
        self.assertEqual(AssociatedLegendre(1, 0, x)[x == 0], 0.0)
        self.assertEqual(AssociatedLegendre(1, 0, x)[x == 4.0], 4.0)

if __name__ == "__main__":
    unittest.main()

