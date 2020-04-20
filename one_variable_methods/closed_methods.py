import pandas as pd
import numpy as np 
from error import relative
from function import Function

pd.set_option("max_colwidth", None)
pd.set_option("display.precision", 16)
pd.set_option('display.max_columns', 500)
pd.set_option("display.width", 2000)

class ClosedMethods:
    '''
    This class contains the followings methods to caltulate 
    some root in an one variable function:
        - Incremental search
        - Bisection
        - False rule
    
    parameters:
        - expr: math expression to work with
        - error: function to calculate the type of error (absolute or relative).
                 relative is set by default.
    '''
    def __init__(self, expr, error = relative):
        self.__function = Function(expr)
        self.error = error


    def incremental_search(self, x0, delta, n):
        '''
        Tries to find out initial values beginning by a random x and going forward 
        incrementally until a sign change happens or n iterations are reached.

        parameters:
            - x0: random x value to begin
            - delta: increment for x in each iteration.
            - n: number of iterations.
        '''
        if n < 1:
            print(f"innapropiate num of iterations = {n}")
        else:
            delta = np.abs(delta) # handling negative values

            fx0 = self.__function.eval(x0)

            data = { "x": np.array([x0], dtype=np.float64),
            "f(x)": np.array([fx0], dtype=np.float64) }

            if fx0 == 0:
                print(f"{x0} is a root")
            else:
                x1 = x0 + delta
                count = 1
                fx1 = self.__function.eval(x1)

                while (fx0 * fx1 > 0) and (count < n):
                    data["x"] = np.append(data["x"], [x1])
                    data["f(x)"] = np.append(data["f(x)"] ,[fx1])

                    x0 = x1 
                    fx0 = fx1
                    x1 += delta
                    fx1 = self.__function.eval(x1)
                    count += 1

                data["x"] = np.append(data["x"], [x1])
                data["f(x)"] = np.append(data["f(x)"] ,[fx1])

                if fx1 == 0:
                    print(f"{x1} is a root")
                elif fx0 * fx1 < 0:
                    print (f" There is a root between {x0} and {x1}")
                else:
                    print(f"Failed in {n} iterations")

            data_frame = pd.DataFrame(data, columns=( "x", "f(x)"), dtype=np.float64)
            data_frame.index.name = "n"
            print(data_frame)
    

    def bisection(self, xi, xs, tolerance, n, false_rule=False):
        '''
        Given a continuos function within [xi, xs] interval and a sign change between it, 
        it tries to find out the closest value to a root by calculating a mean of x.

        parameters:
            - xi: lower end of the interval.
            - xs: higher end of the interval.
            - tolerance: maximun allowed error.
            - n: number of iterations.
            - false_rule: set True if you want to choose false rule method
        '''

        if tolerance < 0 :
            print (f"innapropiate tolerance = {tolerance}")
        elif n < 1:
            print(f"innapropiate num of iterations = {n}")
        else:

            data = {"xi": np.array([],dtype=np.float64), 
            "xs": np.array([],dtype=np.float64),
            "xm" : np.array([],dtype=np.float64), 
            "f(xm)": np.array([],dtype=np.float64), 
            "E": np.array([],dtype=np.float64) }

            fxi = self.__function.eval(xi)
            fxs = self.__function.eval(xs)
            if fxi == 0:
                print (f"{xi} is a root")
            elif fxs == 0:
                print (f"{xs} is a root")
            elif fxi * fxs < 0:
                xm = np.float64(xi - (fxi * (xs - xi)) / (fxs - fxi)) if false_rule else np.mean(np.array((xi, xs), dtype=np.float64))
                fxm = self.__function.eval(xm)
                count = 1
                error = tolerance + 1

                data["xi"] = np.append(data["xi"], [xi])
                data["xs"] = np.append(data["xs"], [xs])
                data["xm"] = np.append(data["xm"], [xm])
                data["f(xm)"] = np.append(data["f(xm)"], [fxm])
                data["E"] = np.append(data["E"], [np.NaN])

                while error > tolerance and fxm != 0 and count < n:
                    if fxi * fxm < 0:
                        xs = xm
                        fxs = fxm
                    else:
                        xi = xm
                        fxi = fxm
                    x_aux = xm 
                    xm = np.float64(xi - (fxi * (xs - xi)) / (fxs - fxi)) if false_rule else np.mean(np.array((xi, xs), dtype=np.float64))
                    fxm = self.__function.eval(xm)
                    error = self.error(xm, x_aux)
                    count += 1

                    data["xi"] = np.append(data["xi"], [xi])
                    data["xs"] = np.append(data["xs"], [xs])
                    data["xm"] = np.append(data["xm"], [xm])
                    data["f(xm)"] = np.append(data["f(xm)"], [fxm])
                    data["E"] = np.append(data["E"], [error])
                if fxm == 0:
                    print (f"{xm} is a root")
                elif error < tolerance:
                    print (f"{xm} is an approximation to a root with tolerance = {tolerance}")
                else:
                    print (f"Failure in {n} iterations")
            else:
                print ("The interval is inappropriate")

            data_frame = pd.DataFrame(data, columns=( "xi", "xs", "xm", "f(xm)", "E"), dtype=np.float64)
            data_frame.index.name = "n"
            print(data_frame)

    def get_function(self):
        return self.__function.get_function()