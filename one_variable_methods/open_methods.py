import pandas as pd
import numpy as np 
from error import relative
from function import Function


pd.set_option("max_colwidth", None)
pd.set_option("display.precision", 16)
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

    def fixed_point(self,expr_g, x0, tolerance, n):
        '''
        It calculates an approximation to a root through finding an intersection
        between y = x and x = g(x) where g(x) is derivated from the function f.

        parameters:
            - x0 : initial point
            - expr_g: function expresion derivated from f
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

            data = {"x": np.array([x0],dtype=np.float64), 
            "f(x)": np.array([fx],dtype=np.float64), 
            "E": np.array([np.NaN],dtype=np.float64) }

            g = Function(expr_g)

            while fx != 0 and error > tolerance and count < n:
                xn = g.eval(x0)
                fx = self.__function.eval(xn)
                error = self.__error(x0,xn)
                x0 = xn
                count += 1

                data["x"] = np.append(data["x"], [xn])
                data["f(x)"] = np.append(data["f(x)"], [fx])
                data["E"] = np.append(data["E"], [error])
            if fx == 0:
                print(f"{x0} is a root")
            
            elif error < tolerance:
                print(f"{x0} is an approximation to a root with tolerance = {tolerance}")
            else:
                print(f"failure in {n} iterations")

            data_frame = pd.DataFrame(data, columns=( "x","f(x)", "E"), dtype=np.float64)
            data_frame.index.name = "n"
            print(data_frame)
            
    def newtons_method(self, x0, tolerance, n):
        '''
        It calculates an approximation to a root through finding an intersection
        between dy_dx and axis x where dy_dx is the first diff of the function f.

            - x0 : initial point
            - tolerance: maximum allowed error
            - n: number of iterations until failure
        parameters:
        '''

        if tolerance < 0:
            print(f"innapropiate tolerance = {tolerance}")
        elif n < 1:
            print(f"innapropiate number of iterations = {n}")
        else:
            fx = self.__function.eval(x0)
            df = self.__function.dfn(x0, 1)
            count = 0
            error = tolerance + 1

            data = {"x": np.array([x0],dtype=np.float64), 
            "f(x)": np.array([fx],dtype=np.float64), 
            "f'(x)": np.array([df],dtype=np.float64), 
            "E": np.array([np.NaN],dtype=np.float64) }

            while fx != 0 and df != 0 and error > tolerance and count < n:
                xn = x0 - (fx / df)
                fx = self.__function.eval(xn)
                df = self.__function.dfn(xn, 1)
                error = self.__error(x0,xn)
                x0 = xn
                count += 1

                data["x"] = np.append(data["x"], [xn])
                data["f(x)"] = np.append(data["f(x)"], [fx])
                data["f'(x)"] = np.append(data["f'(x)"], [df])
                data["E"] = np.append(data["E"], [error])
            if fx == 0:
                print(f"{x0} is a root")
            elif error < tolerance:
                print(f"{x0} is an approximation to a root with tolerance = {tolerance}")
            elif df == 0:
                print(f"{x0} is a possible multiple root")
            else:
                print(f"failure in {n} iterations")

            data_frame = pd.DataFrame(data, columns=( "x", "f(x)", "f'(x)", "E"), dtype=np.float64)
            data_frame.index.name = "n"
            print(data_frame)

    def secant(self, x0, x1, tolerance, n):
        '''
        It calculates an approximation to a root by finding an intersection
        between the  intersection of the secant line (with x0,x1) and axis x where
        x0, x1 are known initial points.

        parameters:
            - x0 : initial point 1
            - x1 : initial point 2
            - tolerance: maximum allowed error
            - n: number of iterations until failure
        '''
        fx0 = self.__function.eval(x0)
        if tolerance < 0:
            print(f"innapropiate tolerance = {tolerance}")
        elif n < 1:
            print(f"innapropiate number of iterations = {n}")
        elif fx0 == 0:
            print (f"{x0} is a root")
        else:
            fx1 = self.__function.eval(x1)
            den = fx1 - fx0
            count = 0
            error = tolerance + 1
            
            data = {"x": np.array([x0, x1],dtype=np.float64), 
            "f(x)": np.array([fx0, fx1],dtype=np.float64), 
            "den": np.array([np.NaN, den],dtype=np.float64), 
            "E": np.array([np.NaN, np.NaN],dtype=np.float64) }

            while fx1 != 0 and error > tolerance and den != 0 and count < n:
                x2 = x1 - (fx1*(x1-x0))/den
                error = self.__error(x1, x2)
                x0 = x1
                fx0 = fx1
                x1 = x2
                fx1 = self.__function.eval(x1)
                den = fx1 - fx0
                count += 1

                data["x"] = np.append(data["x"], [x1])
                data["f(x)"] = np.append(data["f(x)"], [fx1])
                data["den"] = np.append(data["den"], [den])
                data["E"] = np.append(data["E"], [error])
        
            if fx1 == 0:
                print(f"{x1} is a root")
            elif error < tolerance:
                print(f"{x1} is an approximation to a root with tolerance = {tolerance}")
            elif den == 0:
                print(f"{x1} is a possible multiple root")
            else:
                print(f"failure in {n} iterations")

            data_frame = pd.DataFrame(data, columns=( "x", "f(x)", "den", "E"), dtype=np.float64)
            data_frame.index.name = "n"
            print(data_frame)

            
    def multiple_roots(self, x0, tolerance, n):
        '''
        It calculates an approximation to a multiple root by becoming it a simple root
        and applying Newton's method. 

        parameters:
            - x0 : initial point
            - tolerance: maximum allowed error
            - n: number of iterations until failure
        '''

        if tolerance < 0:
            print(f"innapropiate tolerance = {tolerance}")
        elif n < 1:
            print(f"innapropiate number of iterations = {n}")
        else:
            fx = self.__function.eval(x0)
            df1 = self.__function.dfn(x0, 1)
            df2 = self.__function.dfn(x0, 2)
            count = 0
            error = tolerance + 1

            data = {"x": np.array([x0],dtype=np.float64), 
            "f(x)": np.array([fx],dtype=np.float64), 
            "f'(x)": np.array([df1],dtype=np.float64), 
            "f''(x)": np.array([df2],dtype=np.float64), 
            "E": np.array([np.NaN],dtype=np.float64) }

            while fx != 0 and error > tolerance and count < n:
                xn = x0 - (fx*df1) / ( (df1**2) - fx*df2 )
                fx = self.__function.eval(xn)
                df1 = self.__function.dfn(xn, 1)
                df2 = self.__function.dfn(xn, 2)
                error = self.__error(x0,xn)
                x0 = xn
                count += 1

                data["x"] = np.append(data["x"], [xn])
                data["f(x)"] = np.append(data["f(x)"], [fx])
                data["f'(x)"] = np.append(data["f'(x)"], [df1])
                data["f''(x)"] = np.append(data["f''(x)"], [df2])
                data["E"] = np.append(data["E"], [error])
            if fx == 0:
                print(f"{x0} is a root")
            elif error < tolerance:
                print(f"{x0} is an approximation to a root with tolerance = {tolerance}")
            else:
                print(f"failure in {n} iterations")

            data_frame = pd.DataFrame(data, columns=( "x","f(x)","f'(x)","f''(x)","E"), dtype=np.float64)
            data_frame.index.name = "n"
            print(data_frame)

    def get_function(self):
        return self.__function