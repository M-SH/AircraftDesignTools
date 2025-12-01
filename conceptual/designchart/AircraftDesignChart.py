import math

import numpy as np
from matplotlib import pyplot as plt

from conceptual.AircraftPerformance import AircraftPerformance
from conceptual.empiric.Gudmundsson import Gudmundsson
from conceptual.empiric.Raymer import Raymer
from general.Airspeeds import Airspeeds
from general.Convert import Convert
from general.ISAtmosphere import ISAtmosphere
from model.aircraft.Aircraft import Aircraft

if __name__ == "__main__":
    W_S = np.linspace(0.01,30,31)

    cessna162 = Aircraft('Cessna 162 Skycatcher')
    cessna162.ar = 8.0
    cessna162.s_ref = 120 # ft^2
    cessna162.mtom = 1320 # lbf
    cessna162.c_d_min = 0.0333
    cessna162.cl_max = 1.756
    cessna162.cl_to = 0.5
    cessna162.cd_to = 0.038
    cessna162.ground_run_distance = 640 #ft
    cessna162.roc = 14.67 #ft/s
    cessna162.v_cruise = 108 #knots
    cessna162.service_ceiling = 17000 #ft

    e = Raymer.calc_oswald_factor(cessna162.ar)

    # lift-induced drag constant, Gudmundsson (16-6)
    k = 1 / (math.pi * cessna162.ar * e)

    airport_mu = 0.04 # ground friction constant
    airport_rho_imp = Convert().convert(Convert().convert(ISAtmosphere.rho0, Convert.KG, Convert.SLUG), Convert.METER, Convert.FOOT, -3)

    t_w_roll_distance = AircraftPerformance.takeoff_ground_run_distance_imp(W_S, cessna162.ground_run_distance, cessna162.cl_to, cessna162.cl_max, cessna162.cd_to, airport_rho_imp, airport_mu)

    v_upsilon = Gudmundsson.calc_best_roc_speed_single_engine_cas(W_S)
    v_upsilon = Convert().convert(v_upsilon, Convert.KNOT, Convert.FT_PER_S)

    t_w_roc_sealevel = AircraftPerformance.rate_of_climb_imp(W_S, cessna162.roc, cessna162.c_d_min, k, v_upsilon, airport_rho_imp)

    cruise_altitude_ft = 8000
    n = 1.155 # 30 degrees bank angle

    rho = ISAtmosphere().get_density(Convert().convert(cruise_altitude_ft, Convert.FOOT, Convert.METER))
    rho_slugs_per_ft = Convert().convert(Convert().convert(rho, Convert.KG, Convert.SLUG), Convert.METER, Convert.FOOT, -3)
    v_cruise = Convert().convert(cessna162.v_cruise, Convert.KNOT, Convert.FT_PER_S)

    t_w_constant_velocity_turn = AircraftPerformance.constant_velocity_turn_imp(W_S, cessna162.c_d_min, k, n, v_cruise, rho_slugs_per_ft)

    t_w_cruise = AircraftPerformance.cruise_imp(W_S, cessna162.c_d_min, k, v_cruise, rho_slugs_per_ft)

    rho_service_ = ISAtmosphere().get_density(Convert().convert(cessna162.service_ceiling, Convert.FOOT, Convert.METER))
    rho_service_ceiling_imp = Convert().convert(Convert().convert(rho, Convert.KG, Convert.SLUG), Convert.METER, Convert.FOOT, -3)
    v_upsilon = Gudmundsson.calc_best_roc_speed_single_engine_cas(W_S)
    v_upsilon = Convert().convert(v_upsilon, Convert.KNOT, Convert.FT_PER_S)

    # convert to true airspeed at service ceiling altitude
    sigma = ISAtmosphere().get_density_ratio(Convert().convert(cessna162.service_ceiling, Convert.FOOT, Convert.METER))
    v_upsilon_service_ceiling = Airspeeds.to_ktas(Airspeeds.to_keas(v_upsilon)) # TODO

    t_w_service_ceiling = AircraftPerformance.service_ceiling_imp(W_S, cessna162.c_d_min, k, v_upsilon, rho_service_ceiling_imp)

    # plt.plot(WS_ATR_kgm2, TW_ATR_guess, 'o', label='ATR-72 approx point')
    plt.plot(W_S, t_w_roll_distance, label='Ground Roll Distance Take-off ')
    plt.plot(W_S, t_w_roc_sealevel, label='Rate of Climb Take-off')
    plt.plot(W_S, t_w_constant_velocity_turn, label='Constant Velocity Turn Cruise')
    plt.plot(W_S, t_w_cruise, label='Cruise')
    plt.plot(W_S, t_w_service_ceiling, label='Service Ceiling')


    plt.xlabel('Wing loading W/S (kg/m^2)')
    plt.ylabel('Thrust-to-weight ratio T/W')
    plt.title('Constraint Diagram')
    plt.grid(True)
    plt.xlim(0.0, 30.0)
    plt.ylim(0.0, 0.8)
    plt.legend(loc='upper right')

    plt.show()