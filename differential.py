# -------------------------------------------------------------------
# Library Program:
# The following program contains various functions to solve ordinary
# differential equations using various methods.
# -------------------------------------------------------------------

# EULER APPROXIMATE FIRST DIFFERENTIAL METHODS
euler_f = lambda x,h,f: (f(x+h) - f(x)) / h
euler_c = lambda x,h,f: (f(x+h/2) - f(x-h/2)) / h
euler_b = lambda x,h,f: (f(x) - f(x-h)) / h

# EULER APPROXIMATE SECOND DIFFERENTIAL METHODS
euler2_f = lambda x,h,f: (euler_f(x+h,h,f) - euler_f(x,h,f)) / h
euler2_c = lambda x,h,f: (euler_c(x+h/2,h,f) - euler_c(x-h/2,h,f)) / h
euler2_b = lambda x,h,f: (euler_b(x,h,f) - euler_b(x-h,h,f)) / h


# APPROXIMATE SECOND DIFFERENTIAL METHODS
d2_c = lambda x,h,f: (f(x+h) - 2*f(x) + f(x-h)) / h**2