import pandas as pd
import numpy as np 
from error import relative

class ClosedMethods:
    '''
    This class contains the followings methods to caltulate 
    some root in an one variable function:
        - Incremental search
        - Bisection
        - False rule
    
    parameters:
        - function: function to work with
        - error: function to calculate the type of error (absolute or relative).
                 relative is set by default.
    '''
    def __init__(self, function, error = relative):
        self.function = function
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
        fx0 = self.function(x0)

        data = { "x": np.array([x0], dtype=np.float64),
         "f(x)": np.array([fx0], dtype=np.float64) }

        if fx0 == 0:
            print(f"{x0} is root")
        else:
            x1 = x0 + delta
            count = 1
            fx1 = self.function(x1)

            while (fx0 * fx1 > 0) and (count < n):
                data["x"] = np.append(data["x"], [x1])
                data["f(x)"] = np.append(data["f(x)"] ,[fx1])

                x0 = x1 
                fx0 = fx1
                x1 += delta
                fx1 = self.function(x1)
                count += 1

            data["x"] = np.append(data["x"], [x1])
            data["f(x)"] = np.append(data["f(x)"] ,[fx1])

            if fx1 == 0:
                print(f"{x1} is root")
            elif fx0 * fx1 < 0:
                print (f" There is a root between {x0} and {x1}")
            else:
                print(f"Failed in {n} iterations")

        data_frame = pd.DataFrame(data, columns=( "x", "f(x)"), dtype=np.float64)
        data_frame.index.name = "n"
        print(data_frame)
    

    def bisection(self, xi, xs, tolerance, n):
        '''
        Given a continuos function within [xi, xs] interval and a sign change between it, 
        it tries to find out the closest value to a root by calculating a statistical mean of x.

        parameters:
            - xi: lower end of the interval.
            - xs: higher end of the interval.
            - tolerance: maximun allowed error.
            - n: number of iterations.
        '''

        data = {"xi": np.array([],dtype=np.float64), 
         "xs": np.array([],dtype=np.float64),
         "xm" : np.array([],dtype=np.float64), 
         "f(xm)": np.array([],dtype=np.float64), 
         "E": np.array([],dtype=np.float64) }

        fxi = self.function(xi)
        fxs = self.function(xs)
        if fxi == 0:
            print (f"{xi} is root")
        elif fxs == 0:
            print (f"{xs} is root")
        elif fxi * fxs < 0:
            xm = np.mean(np.array((xi, xs), dtype=np.float64))
            fxm = self.function(xm)
            count = 1
            error = tolerance + 1

            data["xi"] = np.append(data["xi"], [xi])
            data["xs"] = np.append(data["xs"], [xs])
            data["xm"] = np.append(data["xm"], [xm])
            data["f(xm)"] = np.append(data["f(xm)"], [fxm])
            data["E"] = np.append(data["E"], [error])

            while error > tolerance and fxm != 0 and count < n:
                if fxi * fxm < 0:
                    xs = xm
                    fxs = fxm
                else:
                    xi = xm
                    fxi = fxm
                x_aux = xm 
                xm = np.mean(np.array((xi, xs), dtype=np.float64))
                fxm = self.function(xm)
                error = self.error(xm, x_aux)
                count += 1

                data["xi"] = np.append(data["xi"], [xi])
                data["xs"] = np.append(data["xs"], [xs])
                data["xm"] = np.append(data["xm"], [xm])
                data["f(xm)"] = np.append(data["f(xm)"], [fxm])
                data["E"] = np.append(data["E"], [error])
            if fxm == 0:
                print (f"{xm} is root")
            elif error < tolerance:
                print (f"{xm} is an approximation to a root with tolerance = {tolerance}")
            else:
                print (f"Failure in {n} iterations")
        else:
            print ("The interval is inappropriate")

        data_frame = pd.DataFrame(data, columns=( "xi", "xs", "xm", "f(xm)", "E"), dtype=np.float64)
        data_frame.index.name = "n"
        print(data_frame)

    def false_rule(self, xi, xs, tolerance, n):
        '''
        Given a continuos function within [xi, xs] interval and a sign change between it, 
        it tries to find out the closest value to a root by calculating a secant line value of x.

        parameters:
            - xi: lower end of the interval.
            - xs: higher end of the interval.
            - tolerance: maximun allowed error.
            - n: number of iterations.
        '''

        data = {"xi": np.array([],dtype=np.float64), 
         "xs": np.array([],dtype=np.float64),
         "xm" : np.array([],dtype=np.float64), 
         "f(xm)": np.array([],dtype=np.float64), 
         "E": np.array([],dtype=np.float64) }

        fxi = self.function(xi)
        fxs = self.function(xs)
        if fxi == 0:
            print (f"{xi} is root")
        elif fxs == 0:
            print (f"{xs} is root")
        elif fxi * fxs < 0:
            xm = xi - (fxi * (xs - xi)) / (fxs - fxi)
            fxm = self.function(xm)
            count = 1
            error = tolerance + 1

            data["xi"] = np.append(data["xi"], [xi])
            data["xs"] = np.append(data["xs"], [xs])
            data["xm"] = np.append(data["xm"], [xm])
            data["f(xm)"] = np.append(data["f(xm)"], [fxm])
            data["E"] = np.append(data["E"], [error])

            while error > tolerance and fxm != 0 and count < n:
                if fxi * fxm < 0:
                    xs = xm
                    fxs = fxm
                else:
                    xi = xm
                    fxi = fxm
                x_aux = xm 
                xm = xi - (fxi * (xs - xi)) / (fxs - fxi)
                fxm = self.function(xm)
                error = self.error(xm, x_aux)
                count += 1

                data["xi"] = np.append(data["xi"], [xi])
                data["xs"] = np.append(data["xs"], [xs])
                data["xm"] = np.append(data["xm"], [xm])
                data["f(xm)"] = np.append(data["f(xm)"], [fxm])
                data["E"] = np.append(data["E"], [error])
            if fxm == 0:
                print (f"{xm} is root")
            elif error < tolerance:
                print (f"{xm} is an approximation to a root with tolerance = {tolerance}")
            else:
                print (f"Failure in {n} iterations")
        else:
            print ("The interval is inappropriate")

        data_frame = pd.DataFrame(data, columns=( "xi", "xs", "xm", "f(xm)", "E"), dtype=np.float64)
        data_frame.index.name = "n"
        print(data_frame)