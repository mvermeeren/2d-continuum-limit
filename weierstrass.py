#!/usr/bin/env python

from sage.symbolic.function import BuiltinFunction
from sage.misc.latex import latex

class Function_p2p(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'wp2p', nargs=1)
        
    def _derivative_(self, z, diff_param):
        return 12*wp(z)*wpp(z)
        
    def _print_latex_(self, z):
        return "\\wp''(%s)" % (latex(z))

wp2p = Function_p2p()

class Function_pp(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'wpp', nargs=1)
        
    def _derivative_(self, z, diff_param):
        return wp2p(z)
        #return 6*wp(z)**2-g2/2
        
    def _print_latex_(self, z):
        return "\\wp'(%s)" % (latex(z))

wpp = Function_pp()

class Function_p(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'wp', nargs=1)
        
    def _derivative_(self, z, diff_param):
        return wpp(z)
        
    def _print_latex_(self, z):
        return "\\wp(%s)" % (latex(z))

wp = Function_p()

class Function_zeta(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'wzeta', nargs=1)
        
    def _derivative_(self, z, diff_param):
        return -wp(z)
        
    def _print_latex_(self, z):
        return "\\zeta(%s)" % (latex(z))

wzeta = Function_zeta()

class Function_sigma(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'wsigma', nargs=1)

    def _derivative_(self, z, diff_param):
        return wzeta(z)*wsigma(z)

    def _print_latex_(self, z):
        return "\\sigma(%s)" % (latex(z))
        
wsigma = Function_sigma()
