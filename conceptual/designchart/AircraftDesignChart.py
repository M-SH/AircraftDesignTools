import math

import numpy as np
from matplotlib import pyplot as plt

from conceptual.AircraftPerformance import AircraftPerformance
from conceptual.empiric.Raymer import Raymer
from general.Convert import Convert
from general.ISAtmosphere import ISAtmosphere
from model.aircraft.Aircraft import Aircraft

if __name__ == "__main__":
    W_S = np.linspace(0,30,31)

    cessna162 = Aircraft('Cessna 162 Skycatcher')
    cessna162.ar = 8.0
    cessna162.s_ref = 120 # ft^2
    cessna162.mtom = 1320 # lbf
    cessna162.c_d_min = 0.333
    cessna162.cl_max = 1.756
    cessna162.cl_to = 0.5
    cessna162.cd_to = 0.038
    cessna162.ground_run_distance = 640 #ft

    e = Raymer.calc_oswald_factor(cessna162.ar)

    # lift-induced drag constant, Gudmundsson (16-6)
    k = 1 / (math.pi * cessna162.ar * e)

    airport_mu = 0.04 # ground friction constant
    airport_rho_imp = Convert().convert(Convert().convert(ISAtmosphere.rho0, Convert.KG, Convert.SLUG), Convert.METER, Convert.FOOT, -3)

    t_w_roll_distance = AircraftPerformance.takeoff_ground_run_distance_imp(10, cessna162.ground_run_distance, cessna162.cl_to, cessna162.cl_max, cessna162.cd_to, airport_rho_imp, airport_mu)

    # plt.plot(WS_ATR_kgm2, TW_ATR_guess, 'o', label='ATR-72 approx point')
    plt.plot(W_S, t_w_roll_distance, label='Take-off Ground Roll Distance')


    plt.xlabel('Wing loading W/S (kg/m^2)')
    plt.ylabel('Thrust-to-weight ratio T/W')
    plt.title('Constraint Diagram')
    plt.grid(True)
    plt.xlim(0.0, 30.0)
    plt.ylim(0.0, 0.8)
    plt.legend(loc='upper right')

    plt.show()