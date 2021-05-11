# -------------------------------------------------------------------
# Library Program:
# The following program contains various functions to solve ordinary
# integrals using various methods.
# -------------------------------------------------------------------

# PRIMITIVE RULES (SINGLE VARIABLE FUNCTIONS)
def rectangular(a,b,N,f): #removed ELs
    '''
    a ; lower-bound for integration
    b ; upper-bound for integration
    N ; number of iterations/subintervals (must be integer)
    f ; function of x (can be user-defined)
    
    The (left/right) rectangular approximation method computes the
    integral of a function by assuming that the area between the
    function and the axis over which the function is being integrated
    can be defined as the area of a series of rectangles, for a
    series of subintervals bound by successive values of x with
    constant step size ;
        xn - xn-1

    For this approximation method, the left/right uppermost corner of
    a rectangle has coordinates computed using the given function,
    then calculates the rectangular area ;
        f(xn) * (xn - xn-1) [for LRAM], or
        f(xn) * (xn+1 - xn) [for RRAM]
    
    These areas are iterated for all subintervals n in N, with the
    sum of the recatangular areas approximating the integral for the
    given function between its lower/upper bounds.
    '''
    dx = (b-a)/N
    Rl = 0
    Rr = 0
    for n in range(0,N):
        Rl += dx * f(a+n*dx)
    for n in range(1,N+1):
        Rr += dx * f(a+n*dx)
    return Rl,Rr

def midpoint(a,b,N,f): #removed ELs
    '''
    a ; lower-bound for integration
    b ; upper-bound for integration
    N ; number of iterations/subintervals (must be integer)
    f ; function of x (can be user-defined)

    The midpoint approximation method operates similarly to the
    rectangular approximation method, in that it approximates the
    area coinciding to the integral of a given function using
    rectangles over a series of subintervals.

    Where the midpoint approximation method differs is that it
    computes the midpoint of two subinterval bounds ;
        Mn = (xn - xn-1) / 2
    
    Once again, like the rectangular method a rectangular area is
    computed using two consecutive subintervals as the x bounds for
    the area and the function of the midpoint as the height ;
        A = f(Mn) * (xn - xn-1)
    
    As for all methods using shapes to approximate the integral of a
    given function within bounds, the areas of the rectangles are
    summed to give a value.
    '''
    dx = (b-a)/N
    M = 0
    for n in range(1,N+1):
        M += f((a+n*dx-a+(n-1)*dx)/2) * dx
    return M

def trapezoid(a,b,N,f): #removed ELs
    '''
    a ; lower-bound for integration
    b ; upper-bound for integration
    N ; number of iterations/subintervals (must be integer)
    f ; function of x (can be user-defined)
    '''
    dx = (b-a)/N
    T = dx/2 * (f(a) + f(b))
    for n in range(1,N):
        T += dx * f(a+n*dx)
    return T

def simpson(a,b,N,f): #removed ELs
    '''
    a ; lower-bound for integration
    b ; upper-bound for integration
    N ; number of iterations/subintervals (must be integer)
    f ; function of x (can be user-defined)
    '''
    dx = (b-a)/N
    S = dx/3 * (f(a) + f(b))
    if N%2 == 0:
        for n in range(1,int(N/2)):
            S += dx/3 * (4*f(a+(2*n-1)*dx) + 2*f(a+2*n*dx))
    if N%2 != 0:
        for n in range(1,int((N+1)/2)):
            S += dx/3 * (4*f(a+(2*n-1)*dx) + 2*f(a+2*n*dx))
    return S
# LARGE UNCERTAINTY AT LOW ITERATIONS (REQUIRES HIGH ITERATIONS TO PRODUCE SIMILAR INTEGRALS TO OTHER METHODS/VERY INCONSISTENT)
def boole(a,b,N,f): #removed ELs
    '''
    a ; lower-bound for integration
    b ; upper-bound for integration
    N ; number of iterations/subintervals (must be integer)
    f ; function of x (can be user-defined)
    '''
    # SECOND TERM (32) IS 2,6,10,...,-10,-6,-2
    # THIRD TERM (12) IS 4,8,12,...,-12,-8,-4
    dx = (b-a)/N
    B = 2*dx/45 * 7*(f(a) + f(b))
    
    Ns_1 = []
    Ns_2 = []
    Ns_3 = []
    
    if N%2 == 0:
        Nrange = int(N/2)
    else:
        Nrange = int((N-1)/2 + 1)
    for n in range(1,Nrange):
        if n%2 == 1:
            Ns_1.append(n)
        if n%4 == 2:
            Ns_2.append(n)
        if n%4 == 0:
            Ns_3.append(n)
    
    for n in Ns_1:
        B += 2*dx/45 * (32*f(a+n*dx) + 32*f(b-n*dx))
    for n in Ns_2:
        B += 2*dx/45 * (12*f(a+n*dx) + 12*f(b-n*dx))
    for n in Ns_3:
        B += 2*dx/45 * (14*f(a+n*dx) + 14*f(b-n*dx))
    return B

def romberg(a,b,f): #removed ELs
    '''
    a ; lower-bound for integration
    b ; upper-bound for integration
    f ; function of x (can be user-defined)
    '''
    test_converge = 1E-9
    dx = (b-a)/2
    CRT = dx/2 * (f(a) + f(b))
    for n in range(1,3):
        CRT += dx * f(a+n*dx)
    R_matrix = [CRT]
    
    i = 2
    while i >= 2:
        dx = (b-a)/(2**(i-1))
        R_j = 0
        for n in range(0,int(2**(i-1)+1)):
            R_j += dx * f(a+n*dx)
        
        Rs = [R_j]
        for j in range(2,i+1):
            R = lambda j: (4**(j-1)*R_j - R_matrix[(j-2)]) / (4**(j-1) - 1)
            Rs.append(R(j))
            R_j = R(j)
            if abs(R(j) - R_j) <= test_converge:
                return Rs[-1]
        R_matrix = Rs
        
        i += 1
# MULTI-VARIABLE FUNCTIONS
def double_rectangular(a,b,c,d,N,M,f): #removed ELs
    '''
    a ; lower-bound for integration over x
    b ; upper-bound for integration over x
    c ; lower-bound for integration over y
    d ; upper-bound for integration over y
    N ; number of x iterations/subintervals (must be integer)
    M ; number of y iterations/subintervals (must be integer)
    f ; function of x and y (can be user-defined)
    '''
    dx = (b-a)/N
    dy = (d-c)/M
    dA = dx*dy
    DR = 0
    for n in range(1,N+1):
        for m in range(1,M+1):
            DR += dA * f(a+n*dx,c+m*dy)
    return DR