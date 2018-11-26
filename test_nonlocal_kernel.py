import unittest
from nonlocal_kernel import Kernel_calculator
import numpy as np

class TestKernelCalculator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.testsample = Kernel_calculator(200,5,8000,8000)

    @classmethod
    def tearDownClass(self):
        pass

    def test_kernel_generator(self):
        second_order_kernel = self.testsample.kernel_generator(2,0.001)
        fourth_order_kernel = self.testsample.kernel_generator(4,0.001)
        np.testing.assert_allclose(np.array(second_order_kernel)[:3], np.array([161.3418,-11.4752,0.8162]),atol=1e-4)
        np.testing.assert_allclose(np.array(fourth_order_kernel)[:3], np.array([175.0393,-17.8644,2.6061]),atol=1e-4)

    def test_frequency_kernel_func(self):
        test_point1 = 24
        test_point2 = 13
        fitting_coeff1 = self.testsample.fitting_kernel_formular(2)[0]
        fitting_coeff2 = self.testsample.fitting_kernel_formular(4)[0]
        
        test_result1 = self.testsample.frequency_kernel_func(test_point1,2,fitting_coeff1)
        test_result2 = self.testsample.frequency_kernel_func(test_point2,2,fitting_coeff1)

        test_result3 = self.testsample.frequency_kernel_func(test_point1,4,fitting_coeff2)
        test_result4 = self.testsample.frequency_kernel_func(test_point2,4,fitting_coeff2)

        self.assertAlmostEqual(test_result1,-12.9050,places=4)
        self.assertAlmostEqual(test_result2,-5.9712,places=4)
        self.assertAlmostEqual(test_result3,-14.1290,places=4)
        self.assertAlmostEqual(test_result4,-6.1044,places=4)



        





if __name__ == '__main__':
    unittest.main()