from math import *


def script(h):
    # ISA CALCULATOR WINAND E4G#
    #############INPUT

    '''
    print("     **** ISA calculator ****")
    print()
    print("1. Calculate ISA for altitude in meters")
    print("2. Calculate ISA for altitude in feet")
    print("3. Calculate ISA for altitude in FL")
    print()
    c = float(input("Enter your choice:"))

    #############CHOICES

    if c == 1.:
        print()
        h = float(input("Enter altitude [m]: "))

    elif c == 2.:
        print()
        h = float(input("Enter altitude [ft]: "))
        h = h * 0.3048

    elif c == 3.:
        print()
        h = float(input("Enter altitude [FL]: "))
        h = h * 30.48

    #############TEMPERATURE

    print()
    print("                        Sea level temperature")
    t0 = (input("Press enter when using standard temperature):"))
    if t0 == "":
        t0 = 288.15
    else:
        t0 = float(t0)
    '''

    t0 = 288.15
    #############CONSTANTS

    h1 = min(h, 11000.)
    h2 = min(h, 20000.)
    h3 = min(h, 32000.)
    h4 = min(h, 47000.)
    h5 = min(h, 51000.)
    h6 = min(h, 71000.)
    h7 = min(h, 84852.)

    g0 = 9.80665
    r = 287.05
    p0 = 101325
    rho0 = 1.225
    a1 = -0.0065
    a2 = 0.0000
    a3 = 0.0010
    a4 = 0.0028
    a5 = 0.0000
    a6 = -0.0028
    a7 = -0.0020

    #############VALUES
    t1 = t0 + h1 * a1
    p1 = p0 * (t1 / t0) ** (-g0 / (a1 * r))
    rho1 = rho0 * (t1 / t0) ** (-g0 / (a1 * r) - 1)

    t2 = t1 + (h2 - h1) * a2
    p2 = p1 * e ** (-g0 / (r * t2) * (h2 - h1))
    rho2 = rho1 * e ** (-g0 / (r * t2) * (h2 - h1))

    t3 = t2 + (h3 - h2) * a3
    p3 = p2 * (t3 / t2) ** (-g0 / (a3 * r))
    rho3 = rho2 * (t3 / t2) ** (-g0 / (a3 * r) - 1)

    t4 = t3 + (h4 - h3) * a4
    p4 = p3 * (t4 / t3) ** (-g0 / (a4 * r))
    rho4 = rho2 * (t4 / t3) ** (-g0 / (a4 * r) - 1)

    t5 = t4 + (h5 - h4) * a5
    p5 = p4 * e ** (-g0 / (r * t5) * (h5 - h4))
    rho5 = rho4 * e ** (-g0 / (r * t5) * (h5 - h4))

    t6 = t5 + (h6 - h5) * a6
    p6 = p5 * (t6 / t5) ** (-g0 / (a6 * r))
    rho6 = rho5 * (t6 / t5) ** (-g0 / (a6 * r) - 1)

    t7 = t6 + (h7 - h6) * a7
    p7 = p6 * (t7 / t6) ** (-g0 / (a7 * r))
    rho7 = rho6 * (t7 / t6) ** (-g0 / (a7 * r) - 1)

    #############

    tt = [t0, t1, t2, t3, t4, t5, t6, t7]
    pp = [p0, p1, p2, p3, p4, p5, p6, p7]
    rr = [rho0, rho1, rho2, rho3, rho4, rho5, rho6, rho7]

    #############sorry

    if h > 84852:
        print()
        print("not yet possible :( ")

    #############CALCULATING

    else:
        if h <= 11000:
            oo = 1
        elif 11000 < h <= 20000:
            oo = 2
        elif 20000 < h <= 32000:
            oo = 3
        elif 32000 < h <= 47000:
            oo = 4
        elif 47000 < h <= 51000:
            oo = 5
        elif 51000 < h <= 71000:
            oo = 6
        elif 71000 < h <= 84852:
            oo = 7
        #print()
        tc = tt[oo] - 273.15
        #print("Temperature : ", round(tt[oo], 2), "K", " (", round(tc, 2), "'C)")
        ps = (pp[oo] / p0) * 100
        #print("Pressure    : ", round(pp[oo], 2), "Pa", "(", round(ps, 2), "%SL)")
        rhos = (rr[oo] / rho0) * 100
        #print("Density     : ", round(rr[oo], 2), "kg/m3", "(", round(rhos, 2), "%SL)")\
    return(rr[oo])