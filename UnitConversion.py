class UC():  # a units conversion class
    def __init__(self):
        """
        This unit converter class is useful for the pipe network and perhaps other problems.
        The strategy is (number in current units)*(conversion factor)=(number desired units), for instance:
            1(ft)*(self.ft_to_m) = 1/3.28084 (m)
            1(in^2)*(self.in2_to_m2) = 1*(1/(12*3.28084))**2 (m^2)
        """

    ft_to_m = 1 / 3.28084
    ft2_to_m2 = ft_to_m ** 2
    ft3_to_m3 = ft_to_m ** 3
    m3_to_ft3 = 1/ft3_to_m3
    ft3_to_L = ft3_to_m3 * 1000
    L_to_ft3 = 1 / ft3_to_L
    in_to_m = ft_to_m / 12
    m_to_in = 1 / in_to_m
    in2_to_m2 = in_to_m ** 2
    m2_to_in2 = 1 / in2_to_m2
    g_SI = 9.80665  # m/s^2
    g_EN = 32.174  # 32.174 ft/s^2
    gc_EN = 32.174  # lbm*ft/lbf*s^2
    gc_SI = 1.0  # kg*m/N*s^2
    lbf_to_kg = 1 / 2.20462
    kg_to_lbf = 1/lbf_to_kg
    lbf_to_N = lbf_to_kg * g_SI
    DeltaF_to_DeltaC = 9.0/5.0
    DeltaC_to_DeltaF = 1/DeltaF_to_DeltaC

    #pressure conversion factors
    pa_to_psi = (1 / (lbf_to_N)) * in2_to_m2
    psi_to_bar = 1.0/14.5038
    bar_to_psi = 1/psi_to_bar

    #energy conversion factors
    kJ_to_btu = 1.0/1.05506
    btu_to_kJ = 1/kJ_to_btu
    kJperkg_to_btuperlb = kJ_to_btu/kg_to_lbf
    btuperlb_to_kJperkg = 1/kJperkg_to_btuperlb

    #entropy conversion factors
    btuperlbF_to_kJperkgC = 4.1868
    kJperkgc_to_btuperlbF = 1/btuperlbF_to_kJperkgC

    #specific volume conversion factors
    ft3perlb_to_m3perkg = ft3_to_m3/lbf_to_kg
    m3perkg_to_ft3perlb = ft3perlb_to_m3perkg

    @classmethod  # this notation allows this method to be directly used from the class by UC.viscosityEnglishToSI
    def viscosityEnglishToSI(cls, mu, toSI=True):
        """
        Converts between lb*s/ft^2 and Pa*s
        :param mu: the viscosity in english units
        :param toSI:  True assumes english in, False assumes SI in
        :return: the viscosity in Pa*s if toSI=True, lb*s/ft^2 if toSI=False
        """
        # (lb*s)/ft^2*((3.3 ft/m)^2)*(1kg/2.2lb)*(9.81m/s^2)->(Pa*s)
        cf = (1 / cls.ft2_to_m2) * (cls.lbf_to_kg) * cls.g_SI
        return mu * cf if toSI else mu / cf

    @classmethod
    def densityEnglishToSI(cls, rho, toSI=True):
        """
        Converts between lb/ft^3 and kg/m^3
        :param rho: specific weight or density
        :param toSI:  True assumes english in, False assumes SI in
        :return: density in SI or EN
        """
        # (lb/ft^3)*((3.3ft/m)^3)*(1kg/2.2lb) -> kg/m^3
        cf = cls.lbf_to_kg / cls.ft3_to_m3
        return rho * cf if toSI else rho / cf

    @classmethod
    def head_to_pressure(cls, h, rho, SI=True):
        """
        Convert from height of column of fluid to pressure in consistent units
        :param h: head in height of fluid (in or m)
        :return: pressure in (psi or Pa)
        """
        if SI:  # p = rho*g*h = g*cf
            cf = rho * cls.g_SI / cls.gc_SI  # kg*m/m^3*s^2
            return h * cf
        else:  # p = rho*g*h = g*cf (h in in)
            cf = rho * cls.g_EN / cls.gc_EN * (1 / 12) ** 2  # (lbm*ft/ft^3*s^2)(lbf*s^2/lbm*ft)(ft^2/in^2)
            return h * cf
        # convert m of water to psi
        # (m)*(3.3*12in/m)*rho(kg/m^3)*(2.2lb/kg)*(1m/(3.3*12in))^3
        psi = p * cls.rho * 2.2 / ((3.3 * 12) ** 2)
        return psi

    @classmethod
    def m_to_psi(cls, h, rho):
        """
        For converting from height of fluid to psi
        :param h: height of fluid in m
        :param rho: density of fluid in kg/m^3
        :return: pressure in psi
        """
        return cls.head_to_pressure(h, rho) * cls.pa_to_psi

    @classmethod
    def psi_to_m(cls, p, rho):
        """
        For converting from psi to height of fluid.
        first convert psi to pa
        :param p: pressure in psi
        :param rho: density of fluid in kg/m^3
        :return: height of fluid in m
        """
        pa = p / cls.pa_to_psi
        h = pa / (rho * cls.g_SI)
        return h

    @classmethod
    def C_to_F(cls, T):
        return T*9.0/5.0+32

    @classmethod
    def F_to_C(cls, T):
        return (T-32)*5.0/9.0