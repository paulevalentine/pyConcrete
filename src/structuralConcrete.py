import math
import numpy as np
import sympy as sp
from IPython.core.display_functions import display


class Concrete(object):
    """ Properties of structural concrete design """

    def __init__(self, fck):
        self.fck = fck # cylinder compressive strength of concrete
        self.unitWeight = 25 # unit weight of concrete in kN/m3

        # set parameters as derived from table 3.1 of BS EN 1992-1-1
        self.fcm = self.calc_mean_compressive_strength()

        # set the value of the mean tensile strength fctm
        self.fctm = self.calc_mean_tensile_strength()

        # Young's modulus in MPa
        #self.Ecm = 22 * (self.fcm / 10)**0.30 *10**3
        self.Ecm = self.calc_elastic_modulus()

    # --- standard calculations --- #
    def calc_mean_compressive_strength(self)->float:
        fcm, fck = sp.symbols('f_cm, f_ck')
        fcm_eq = sp.Eq(fcm, fck + 8)
        fcm_val = fcm_eq.subs({fck:self.fck}).evalf(3)
        display(fcm_eq, fcm_val)
        return fcm_val.rhs

    def calc_mean_tensile_strength(self)->float:
        fctm, fck = sp.symbols('f_ctm, f_ck')
        if self.fck <= 50:
            fctm_eq = sp.Eq(fctm, 0.30 * fck**sp.Rational(2,3))
            fctm_val = fctm_eq.subs({fck: self.fck}).evalf(3)
            display(fctm_eq, fctm_val)
            return fctm_val.rhs
        else:
            fcm = sp.symbols('f_cm')
            fctm_eq = sp.Eq(fctm, 2.12 * sp.log(1+fcm/10))
            fctm_val = fctm_eq.subs({fcm:self.fcm}).evalf(3)
            display(fctm_eq, fctm_val)
            return fctm_val.rhs

    def calc_elastic_modulus(self)->float:
        Ecm, fcm = sp.symbols('E_cm, f_cm')
        Ecm_eq = sp.Eq(Ecm, 22*(fcm/10)**0.30 *10**3)
        Ecm_val = Ecm_eq.subs({fcm: self.fcm}).evalf(3)
        display(Ecm_eq, Ecm_val)
        return Ecm_val.rhs

    # --- Time dependent tensile and compressive strength properties --- #
    def bcc(self,s : float, t : float)-> float:
        """ Time function for concrete strength """
        return math.e**(s * (1 - (28/t)**0.50))

    def fcmt (self, s : float, t : float)->float:
        """ Mean compressive strength as a function of time """
        return self.fcm * self.bcc(s,t)

    def fckt(self, s:float, t)->np.ndarray:
        """ Characteristic compressive strength as a function of time """
        # Convert to numpy arrays to handle conditional logic across elements
        t = np.asanyarray(t)

        # Use np.select for multiple conditions
        condlist = [t < 3, (t >= 3) & (t < 28)]
        choicelist = [0, self.fcmt(s, t) - 8]

        return np.select(condlist, choicelist, default=self.fck)