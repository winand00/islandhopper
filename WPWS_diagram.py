from math import *
import numpy as np


def dragcoef(C_D_0, C_L, A, e):
    C_D = C_D_0 + (C_L ** 2) / (pi * A * e)
    return C_D


def Stallload(rho, C_L_max, V_s):
    WS = 0.5 * rho * V_s ** 2 * C_L_max
    return WS


def TOPcalc(S_to):
    a = 0.0577
    b = 8.6726
    c = 0
    TOP = (-b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)


def cruise_perf(x, sigma, n_p, power_setting, cruise_fraction, C_D_0, W, S, rho, V, A, e):
    y = power_setting / cruise_fraction * n_p * sigma ** 0.75 * (
                C_D_0 * 0.5 * rho * V ** 3 / (cruise_fraction * W / S) + (cruise_fraction * W / S) * 1 / (
                    pi * A * e * 0.5 * rho * V)) ** (-1)
    return (x, y)


def climb_rate(x, n_p, c, A, e, C_D_0, rho):
    y = n_p / (c + ((sqrt(x) * sqrt(2 / rho)) / (1.345 * (A * e) ** 0.75 / C_D_0 ** 0.25)))
    return (x, y)


def climb_gradient(x, n_p, cV, C_DC_L, C_L, rho):
    y = n_p / (sqrt(x) * (cV + C_DC_L) * sqrt(2 / rho / C_L)
    return (x, y)
