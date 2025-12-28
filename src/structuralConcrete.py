import math
import numpy as np

class Concrete(object):
    """ Properties of structural concrete design """

    def __init__(self, fck):
        self.fck = fck # cylinder compressive strength of concrete
        self.unitWeight = 25 # unit weight of concrete in kN/m3

        # set parameters as derived from table 3.1 of BS EN 1992-1-1
        self.fcm = fck + 8  # in MPa

        # set the value of the mean tenile strength fctm
        if self.fck <= 50:
            self.fctm = 0.30 * self.fck ** (2 / 3)
        else:
            self.fctm = 2.12 * math.log(1 + self.fcm / 10)

        # set the characteristic value for the tensile strength
        self.fctk = 0.70 * self.fctm

        # Young's modulus in MPa
        self.Ecm = 22 * (self.fcm / 10)**0.30 *10**3

    # Time dependent tensile and compressive strength
    def bcc(self,s : float, t : float)-> float:
        """ Time function for concrete strength """
        return math.e**(s * (1 - (28/t)**0.50))

    def fcmt (self, s : float, t : float)->float:
        """ Mean compressive strength as a function of time """
        return self.fcm * self.bcc(s,t)

    def fckt(self, s:float, t:float)->np.ndarray:
        """ Characteristic compressive strength as a function of time """
        # Convert to numpy arrays to handle conditional logic across elements
        t = np.asanyarray(t)

        # Use np.select for multiple conditions
        condlist = [t < 3, (t >= 3) & (t < 28)]
        choicelist = [0, self.fcmt(s, t) - 8]

        return np.select(condlist, choicelist, default=self.fck)
