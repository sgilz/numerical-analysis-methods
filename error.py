from numpy import abs as ABS

def absolute(xi, xf):
    '''
    Calculates absolute error
    params:
        - xi: previous value
        - xf: current value
    '''

    if xi == None or xf == None:
        raise AttributeError("values xi and xf cannot be None") 

    return ABS(xi - xf)

def relative(xi, xf):
    '''
    Calculates relative error
    params:
        - xi: previous value
        - xf: current value

    warning:
    If xf (current value) is zero, a ZeroDivisionError is raised
    '''

    if xf == 0:
        raise ZeroDivisionError("Current value (xf) cannot be zero")
    if xi == None or xf == None:
        raise AttributeError("values xi and xf cannot be None")
    
    return ABS( (xi - xf)/ xf )