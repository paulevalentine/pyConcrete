import structuralConcrete
import math
import matplotlib.pyplot as plt
import numpy as np

class RectBeam(structuralConcrete.Concrete):
    """ Reinforced concrete beam calculations """
    def __init__(self, fck, b, h, cover, link_diam, bar_diam):
        super().__init__(fck)
        self.b = b
        self.h = h
        self.cover = cover
        self.link_diameter = link_diam
        self.bar_diameter = bar_diam
        self.effective_depth = self.h - self.cover - self.link_diameter - self.bar_diameter / 2
        self.gc = 1.5 # partial factor on material strength
        self.fyk = 500 # characteristic strength for high yield steel
        self.acc = 0.85
        delta = 1
        kMax = 0.6 * delta - 0.18 * delta ** 2 - 0.21
        self.m_max = kMax * (self.b * self.effective_depth ** 2 * self.fck) * 10 ** -6
    def bending_steel_required(self, moment: float, p: str = "yes")->float:
        """ This function does not establish the steel if compression steel is required """

        if moment > self.m_max:
            print(
                f'Compression steel is require. The maximum moment the section can take without compression reinforcement is {m_max :.2f}kNm')
            return 0

        k = moment * 10 ** 6 / (self.b * self.effective_depth ** 2 * self.fck)

        # calculate the level arm
        z = min(self.effective_depth / 2 * (1 + (1 - 3.53 * k) ** 0.50), 0.95 * self.effective_depth)

        # calculate the requited bending steel
        ast = moment * 10 ** 6 / (0.87 * 500 * z)
        if p == "yes":
            print(f"The area of longitudinal bending steel required = {ast:.2f}mm2")
        return ast

    def shear_cap_no_links(self, ast_prov : float)->float:
        """ Check the shear capacity of a section if no shear reinforcement is provided """
        p1 = min(ast_prov / (self.b * self.effective_depth), 0.02)
        kv = min(1 + (200 / self.effective_depth) ** 0.50, 2)
        crdc = 0.18/self.gc
        vrdc = max(crdc*kv*(100*p1*self.fck)**(1/3),
                   0.035*kv**1.5*self.fck**0.50)

        vcap = self.b * self.effective_depth * vrdc / 1000
        print(f"The shear capcity of the section with no links = {vcap:.2f}kN")

        return vcap

    def shear_cap_with_links(self, asw : float, sv : float, theta : float)->float:
        """ check the shear capacity of a section if shear reinforcement is provided """
        # check to make sure max link spacing provided
        if(sv > 0.75*self.effective_depth):
            print("The links spacing specified is greater than the max spacing")
        f = theta * math.pi/180

        v = asw * self.effective_depth * 0.90 * 0.87 * self.fyk * (1 / math.tan(f)) / (sv * 1000)

        fcd = self.acc*self.fck/self.gc

        self.vMax = self.b * self.effective_depth * 0.90 * 0.60 * (1 - self.fck / 250) * fcd / (1 / math.tan(f) + math.tan(f)) * (1 / 1000)

        vcap = min(v, self.vMax)
        print(f"The shear capacity of the section with links = {vcap:.2f} kN")

        return vcap

    def span_to_depth(self, ast_req:float, ast_prov:float, comp_steel:float, K:float, span:float)->float:
        """ method to calculate the allowable span to depth ratio """
        p0 = self.fck** 0.50 * 10 ** -3
        p = ast_req / (self.b * self.effective_depth)
        pc = comp_steel / (self.b * self.effective_depth)

        if p < p0:
            sdr = K * ((11 + 1.5 * self.fck**0.5 * (p0 / p) + 3.2 *
                        (p0 / p - 1) ** 1.5) * (500 / (500 * ast_req / ast_prov)))
        else:
            sdr = K * ((11 + 1.5 * self.fck ** 0.5 * (p0 / (p - pc)) + (1 / 12) * self.fck
                        ** 0.50 * (pc / p0) ** 0.50) * (500 / (500 * ast_req / ast_prov)))

        print(f"Required maximum span to depth ratio = {sdr:.2f}")
        ld = span * 10 ** 3 / self.effective_depth
        print(f"Actual spand to depth ratio = {ld:.2f}")
        return sdr

    def max_torsion_cap(self, theta:float, asw:float, sv:float, ast_prov:float, Ved:float, Ted:float)->None:
        vCap = self.shear_cap_with_links(asw, sv, theta)  # call to assign vMax
        t = (self.b * self.h) / (2 * (self.b + self.h))
        Ak = (self.h - t) * (self.b - t)
        uk = 2 * ((self.h - t) + (self.b - t))
        tefi = max(t, 2 * (self.h - self.effective_depth))
        fcd = self.acc * self.fck / self.gc
        acw = 1
        v = 0.6 * (1 - (self.fck / 250))
        th = theta * math.pi / 180  # convert to radians

        #  check the limit from the compressive struts
        TRdmax = 2 * v * acw * fcd * Ak * tefi * math.sin(th) * math.cos(th) * 10 ** -6
        print(f"Max torsional capacity of the section = {TRdmax:.2f} kNm")
        print(f"The maximum shear capacity of the section = {self.vMax:.2f} kN")
        uRatioA = (Ted / TRdmax) + (Ved / self.vMax)
        print(f"Concrete strut unity ratio = {uRatioA:.2f}")

        # calculate the nominal torsional capacity of the section
        fctd = self.fctk / self.gc
        TRdc = 2 * Ak * fctd * t * 10 ** -6
        print(f"Torsional cracking moment = {TRdc:.2f} kNm")
        v_cap_no_links = self.shear_cap_no_links(ast_prov)
        print(f"Unreinforced shear capacity of section = {v_cap_no_links:.2f} kN")
        uRatioB = (Ted / TRdc) + (Ved / v_cap_no_links)
        print(f"Minimum reinforcement unity check = {uRatioB:.2f}")
        if uRatioB > 1:
            print(f"Design torsion reinforcement")
            ti = Ted * 10 ** 6 / (2 * Ak * t)
            print(f"Shear stress in effective wall = {ti:.2f} MPa")
            Vt = ti * t * self.effective_depth * 10 ** -3
            VTotal = Ved / 2 + Vt
            print(f"Shear force in a vertical wall = {Vt:.2f} kN")
            print(f"Total shear force in the vertical wall = {VTotal:.2f} kN")
            vWallCap = self.shear_cap_with_links(asw / 2, sv, theta)
            print(f"Shear capacity of the wall = {vWallCap:.2f} kN")
            if vWallCap > VTotal:
                print("Torsional capacity okay")

                # calculate the additional longtinudinal reinforcement required
                # in each effective wall
                Ftd = 0.50 * VTotal * (1 / math.tan(th))
                print(f"Additonal tension in longitudinal reinforcement = {Ftd:.2f} kN")
                astAdd = Ftd * 10 ** 3 / (0.87 * self.fyk)
                print(f"Additional longitudinal reinforcement required = {astAdd:.2f} mm2")
            else:
                print('FAIL')
        else:
            print('No designed shear reinforcement required')

    # ------ print capacities ----- #

    def print_bending_capacity_curve(self)->None:
       y = np.linspace(0,self.m_max,100)
       x = np.array([])
       for vals in y:
           ast = self.bending_steel_required(vals, "n")
           x = np.append(x,ast)
       plt.plot(x,y, label="Capacity curve")
       plt.xlabel("Area of bending steel (mm2)")
       plt.ylabel("ULS Bending Capacity (kNm)")
       a = math.pi * 16**2
       b = math.pi * 20**2
       c = math.pi * 25**2
       d = math.pi * 32**2
       plt.axvline(a, color='r', linestyle='--', label="4H16")
       plt.axvline(b, color='g', linestyle='--', label="4H20")
       plt.axvline(c, color='b', linestyle='--',label="4H25")
       plt.axvline(d, color='y', linestyle='--', label="4H32")
       plt.title(f"{self.b}mmx{self.h}mm Beam in Bending")
       plt.grid()
       plt.legend()
       plt.show()

    def print_shear_capacity_curve(self)->None:
        print("This does nothing yet")
        #todo write code for the shear curve




