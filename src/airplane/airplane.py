from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

import warnings

from ..kernel import gauss_chebyshev as gc


"""
            COSAS POR HACER
        - Revisar las unidades para acabar con newtons 
        - Impedir que se usen lifts negativos exagerados
        - Si da tiempo, añadir las series de Fourier
        - En vez de fourier, se pueden hacer varias aproximaciones polinomicas. Hacerlo todo en una funcion que mida el length de a
        - Si da tiempo, penalizar el drag inducido para evitar muchisimo lift en la raiz o el lift en la punta
"""


class Airplane:
    def __init__(self, mass: float, wingspan: float):

        self.mass = mass
        self.weight = mass * 9.8
        self.wingspan = wingspan    



    # NOTE NOTE    Revisar unidades    NOTE NOTE

    def lift_poli(self, a: list[float]) -> float:
        """
        Integrate the lift distribution in the wing. The moment distribution can be a polynomial of variable degree.

        Args:
        a: list of floats. Coefficients of the polynomial in the form a[i]*x**(i+1)

        Returns:
        float: the result of the integration
        """
        return self._integrate(self._polinomial, -1, 1, a)
    


    def lift_poli_fortran(self, n : int, a: list[float]) -> float:
        return self._fortran_integrate_lift(n, a)
    


    def moment_poli(self, a: list[float]) -> float:
        """
        Integrate the moment distribution in the wing. The moment distribution can be a polynomial of variable degree.

        Args:
        a: list of floats. Coefficients of the polynomial in the form a[i]*x**(i+1)

        Returns:
        float: the result of the integration
        """
        return self._integrate(self._polinomial_x, -1, 1, a) * self.wingspan/2
    

    def moment_poli_fortran(self, n : int, a: list[float]) -> float:
        """
        Integrate the moment distribution in the wing. The moment distribution is a 3th degree polinomial

        Args:
        a: list of 3 floats. Coefficients of the polinomial in the form a[0]*x**3 + a[1]*x**2 + a[2]*x

        Returns:
        float: the result of the integration
        """
        return self._fortran_integrate_moment(n, a) * self.wingspan/2
    

    def lift_fourier(self, a: list[float]) -> float:
        raise NotImplementedError("Fourier series not implemented yet")

    def moment_fourier(self, a: list[float]) -> float:
        raise NotImplementedError("Fourier series not implemented yet")
    

    def abs_lift(self, a: list[float]) -> float:
        return self._integrate(self._abs_polinomial, -1, 1, a)

    

    def plot_lift(self, a: list[float]):
        x = np.linspace(-1, 1, 100)
        y = self._polinomial(x, a)
        plt.plot(x, y)
        plt.title("Distribucion de sustentacion")
        plt.xlabel("Posicion en la envergadura")
        plt.ylabel("Sustentacion [N]")
        plt.show()

    def plot_moment(self, a: list[float]):
        x = np.linspace(-1, 1, 100)
        y = self._polinomial_x(x, a)
        plt.plot(x, y)
        plt.title("Distribucion de momentos")
        plt.xlabel("Posición en la envergadura")
        plt.ylabel("Momento [N*m]")
        plt.show()


    def _polinomial(self, x : float | np.ndarray, a: list[float]) -> float:
        """
        Polinomial. Load distribution in the wing.

        Args:
        x: float or np.ndarray in the interval [-1, 1]
        a: list of 4 floats. Coefficients of the polinomial in the form a[0]*x + a[1]*x**2 + ...
        """

        if not self._check_x(x):
            raise ValueError("x values must be in the interval [-1, 1]")
        
        length = len(a)
        result = 0

        for i in range(length):
            result += a[i]*x**(i+1)

        result *= np.sqrt(1-x**2)

        return result
    

    def _polinomial_x(self, x : float, a: list[float]) -> float:
        """
        Moment of the polinomial.
        """
        return self._polinomial(x, a)*(x+1) # NOTE: x+1 se pone porque la integral va de -1 a 1. Revisar si es correcto
    

    def _abs_polinomial(self, x : float | np.ndarray, a: list[float]) -> float:
        return abs(self._polinomial(x, a))
            

        

    def _integrate(self, f : any, min : float, max : float, coefs: list[float]) -> float:
        """
        Integrate a function with the given coefficients. It is designed to 
        integrate private functions of the class.

        Args:
        f: function to integrate
        min: minimum value of the interval
        max: maximum value of the interval
        coefs: list of coefficients of the function. Example: [a, b, c] for f(x) = a*x^2 + b*x + c

        Returns:
        float: the result of the integration
        """
        return integrate.quad(f, min, max, args=(coefs))[0]
    
    def _fortran_integrate_lift(self, n : int, a: list[float]) -> float:
        return gc.run_lift(n, a)
    
    def _fortran_integrate_moment(self, n : int, a: list[float]) -> float: 
        return gc.run_mom(n, a)


    def _check_x(self, x: float) -> bool:
        """
        Check if x is in the interval [0, 1]
        """
        if isinstance(x, (float, int)):
            return x<= 1 and x>= -1
        elif isinstance(x, np.ndarray):
            return (x<= 1).all() and (x>= -1).all()
        else:
            raise TypeError(f"x must be a float, int or a numpy array, but got {type(x)}")
        
    
    def _check_length(self, a: list[float], length: int) -> bool:
        """
        Check if the length of the list is equal to the length given
        """
        return len(a) == length
    



    ############### DEPRECATED METHODS ####################

    def _4polinomial(self, x : float | np.ndarray, a: list[float]) -> float:
        """
        4th degree polinomial. Load distribution in the wing.

        Args:
        x: float or np.ndarray in the interval [0, 1]
        a: list of 4 floats. Coefficients of the polinomial in the form a[0]*x**4 + a[1]*x**3 + a[2]*x**2 + a[3]*x
        """
        warnings.warn("This method is deprecated. Use _polinomial instead", DeprecationWarning)

        if not self._check_x(x):
            raise ValueError("x values must be in the interval [0, 1]")
        if not self._check_length(a, 4):
            raise ValueError(f"a must have length 4, but got {len(a)}")
        
        return (a[0]*x**4 + a[1] * x**3 + a[2] * x**2 + a[3] * x)*np.sqrt(1-x**2)
    

    def _4polinomial_x(self, x : float, a: list[float]) -> float:
        """
        Moment of the 4th degree polinomial.
        """
        warnings.warn("This method is deprecated. Use _polinomial_x instead", DeprecationWarning)
        return self._4polinomial(x, a)*(x+1) # NOTE: x+1 se pone porque la integral va de -1 a 1. Revisar si es correcto


        
