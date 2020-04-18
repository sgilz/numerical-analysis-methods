from sympy import Symbol, diff
from sympy.parsing.sympy_parser import parse_expr
import numpy as np

class Function:
    def __init__(self, expression, var='x'):
        self.__var = Symbol(var)
        self.__function = parse_expr(expression)

    def eval(self, x):
        return np.float64(self.__function.subs(self.__var, x).evalf())

    def dfn(self, x, n):
        if n < 1:
            raise AttributeError(f"{n} is an innapropiate number of iterations")
        df = self.__function
        var = self.__var
        for _ in range(n):
            df = diff(df, var)
        return np.float64(df.subs(var, x).evalf())

    def get_function(self):
        return self.__function