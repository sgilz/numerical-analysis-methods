import pandas as pd
import numpy as np 
from error import relative
from function import Function


pd.set_option("display.max_colwidth", -1)
pd.set_option('display.max_columns', 500)
pd.set_option("display.width", 2000)

class OpenMethods:
    '''
    This class contains the followings methods to caltulate 
    some root in an one variable function:
        - Fixed point
        - Newton's method
        - Secant
        - Multiple roots
    
    Args:
        - function f: function to work with
        - error: function to calculate the type of error (absolute or relative).
                 relative is set by default.
    '''
    
    def __init__( self, expression, error = relative):
        self.__function = Function(expression)
        self.__error = error

    def fixed_point(self, x0, g, tolerance, n):
        '''
        It calculates an approximation to a root through find an intersection
        between y = x and x = g(x) where g(x) is derivated from the function f.

        parameters:
            - x0 : initial point
            - function g: function derivated from f
            - tolerance: maximum allowed error
            - n: number of iterations until failure
        '''
        
        if tolerance < 0:
            print(f"innapropiate tolerance = {tolerance}")
        elif n < 1:
            print(f"innapropiate number of iterations = {n}")
        else:
            fx = self.__function.eval(x0)
            count = 0
            error = tolerance + 1
            while fx != 0 and error > tolerance and count < n:
                xn = g(x0)
                fx = self.__function.eval(xn)
                error = self.__error(x0,xn)
                x0 = xn
                count += 1
            if fx == 0:
                print(f"{x0} is a root")
            
            elif error < tolerance:
                print(f"{x0} is an approximation to a root with tolerance = {tolerance}")
            else:
                print(f"failure in {n} iterations")
            


    def newtons_method(self):
        pass
    def secant(self):
        pass
    def multiple_roots(self):
        pass