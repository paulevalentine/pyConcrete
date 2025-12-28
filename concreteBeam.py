import structuralConcrete
import math

class Beam(structuralConcrete.Concrete):
    """ Reinforced concrete beam calculations """
    def __init__(self, fck, b, h, cover, linkDiam, barDiam):
        super().__init__(fck)
        self.b = b
        self.h = h
        self.cover = cover
        self.linkDiam = linkDiam
        self.barDiam = barDiam
        self.effectiveDepth = self.h - self.cover - self.linkDiam - self.barDiam / 2

        self.gc = 1.5 # partial factor on material strength
        self.fyk = 500 # characteristic strength for high yield steel
        self.acc = 0.85

    def bendingSteelReq(self, moment: float)->float:
        """ This function does not establish the steel if compression steel is required """

        delta = 1
        kMax = 0.6 * delta - 0.18 * delta ** 2 - 0.21
        m_max = kMax * (self.b * self.effectiveDepth ** 2 * self.fck) * 10 ** -6
        if moment > m_max:
            print(
                f'Compression steel is require. The maximum moment the section can take without compression reinforcement is {m_max :.2f}kNm')
            return 0

        k = moment * 10 ** 6 / (self.b * self.effectiveDepth ** 2 * self.fck)

        # calculate the level arm
        z = min(self.effectiveDepth / 2 * (1 + (1 - 3.53 * k) ** 0.50), 0.95 * self.effectiveDepth)

        # calculate the requited bending steel
        ast = moment * 10 ** 6 / (0.87 * 500 * z)
        print(f"The area of longitudinal bending steel required = {ast:.2f}mm2")

        return ast

    def shearCapNoLinks(self, astProv : float)->float:
        """ Check the shear capacity of a section if no shear reinforcement is provided """
        p1 = min(astProv/(self.b*self.effectiveDepth), 0.02)
        kv = min(1+(200/self.effectiveDepth)**0.50, 2)
        cRdc = 0.18/self.gc
        vRdc = max(cRdc*kv*(100*p1*self.fck)**(1/3),
                   0.035*kv**1.5*self.fck**0.50)

        VCap = self.b*self.effectiveDepth * vRdc/1000
        print(f"The shear capcity of the section with no links = {VCap:.2f}kN")

        return VCap

    def shearCapLinks(self, asw : float, sv : float, theta : float)->float:
        """ check the shear capacity of a section if shear reinforcement is provided """
        # check to make sure max link spacing provided
        if(sv > 0.75*self.effectiveDepth):
            print("The links spacing specified is greater than the max spacing")
        f = theta * math.pi/180

        v = asw*self.effectiveDepth*0.90*0.87*self.fyk*(1/math.tan(f))/(sv*1000)

        fcd = self.acc*self.fck/self.gc

        self.vMax = self.b*self.effectiveDepth*0.90*0.60 * (1-self.fck/250)*fcd/(1/math.tan(f)+math.tan(f))*(1/1000)

        vCap = min(v, self.vMax)
        print(f"The shear capacity of the section with links = {vCap:.2f} kN")

        return vCap

    def spanToDepth(self, astReq:float, astProv:float, compSteel:float, K:float, span:float)->float:
        """ method to calculate the allowable span to depth ratio """
        p0 = self.fck** 0.50 * 10 ** -3
        p = astReq / (self.b * self.effectiveDepth)
        pc = compSteel / (self.b * self.effectiveDepth)

        if p < p0:
            sdr = K * ((11 + 1.5 * self.fck**0.5 * (p0 / p) + 3.2 *
                        (p0 / p - 1) ** 1.5) * (500 / (500 * astReq / astProv)))
        else:
            sdr = K * ((11 + 1.5 * self.fck ** 0.5 * (p0 / (p - pc)) + (1 / 12) * self.fck
                        ** 0.50 * (pc / p0) ** 0.50) * (500 / (500 * astReq / astProv)))

        print(f"Required maximum span to depth ratio = {sdr:.2f}")
        ld = span * 10 ** 3 / self.effectiveDepth
        print(f"Actual spand to depth ratio = {ld:.2f}")
        return sdr

    def maxTorsionCap(self, theta:float, asw:float, sv:float, astProv:float, Ved:float, Ted:float)->None:
        vCap = self.shearCapLinks(asw, sv, theta)  # call to assign vMax
        t = (self.b * self.h) / (2 * (self.b + self.h))
        Ak = (self.h - t) * (self.b - t)
        uk = 2 * ((self.h - t) + (self.b - t))
        tefi = max(t, 2 * (self.h - self.effectiveDepth))
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
        vCapNoLinks = self.shearCapNoLinks(astProv)
        print(f"Unreinforced shear capacity of section = {vCapNoLinks:.2f} kN")
        uRatioB = (Ted / TRdc) + (Ved / vCapNoLinks)
        print(f"Minimum reinforcement unity check = {uRatioB:.2f}")
        if uRatioB > 1:
            print(f"Design torsion reinforcement")
            ti = Ted * 10 ** 6 / (2 * Ak * t)
            print(f"Shear stress in effective wall = {ti:.2f} MPa")
            Vt = ti * t * self.effectiveDepth * 10 ** -3
            VTotal = Ved / 2 + Vt
            print(f"Shear force in a vertical wall = {Vt:.2f} kN")
            print(f"Total shear force in the vertical wall = {VTotal:.2f} kN")
            vWallCap = self.shearCapLinks(asw / 2, sv, theta)
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
