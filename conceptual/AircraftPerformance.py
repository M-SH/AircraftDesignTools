import numpy as np
from numpy import ndarray

from general.Airspeeds import Airspeeds
from general.Convert import Convert
from general.ISAtmosphere import ISAtmosphere
from general.PhysicalConstant import PhysicalConstant
from general.Physics import Physics


class AircraftPerformance:
    @staticmethod
    def takeoff_ground_run_distance(wing_loading: float | ndarray, ground_run_distance: float, cl_to, cl_max, cd_to: float, rho: float = 1.225, my: float = 0.04) -> float | ndarray:
        term1 = (1.21 / (PhysicalConstant.G * rho * cl_max * ground_run_distance)) ** wing_loading
        term2 = (0.605 / cl_max) * (cd_to - my * cl_to)
        term3 = my

        thrust_to_weight = term1 + term2 + term3

        return thrust_to_weight

    @staticmethod
    def takeoff_ground_run_distance_imp(wing_loading: float | ndarray, ground_run_distance: float, cl_to, cl_max, cd_to: float, rho: float = 1.225, my: float = 0.04) -> float | ndarray:
        g_imp = Convert().convert(PhysicalConstant.G, Convert.METER, Convert.FOOT)

        term1 = (1.21 * wing_loading / (g_imp * rho * cl_max * ground_run_distance))
        term2 = (0.605 / cl_max) * (cd_to - my * cl_to)
        term3 = my

        thrust_to_weight = term1 + term2 + term3

        return thrust_to_weight

    @staticmethod
    def rate_of_climb_imp(wing_loading: float | ndarray, roc: float, c_d_min: float, k: float, v_upsilon_kcas: float, altitude_m: float) -> float | ndarray:
        rho = ISAtmosphere().get_density(altitude_m)
        rho_imp = Convert().convert(Convert().convert(rho, Convert.KG, Convert.SLUG), Convert.METER, Convert.FOOT, -3)

        v_upsilon_ktas = Airspeeds.to_ktas(Airspeeds.to_keas(v_upsilon_kcas, altitude_m), altitude_m)
        v_upsilon_ft_per_s = Convert().convert(v_upsilon_ktas, Convert.KNOT, Convert.FT_PER_S)
        q = Physics.calc_dynamic_pressure(rho_imp, v_upsilon_ft_per_s)

        term1 = roc / v_upsilon_ft_per_s
        term2 = (q / wing_loading) * c_d_min
        term3 = (k / q) * wing_loading

        thrust_to_weight = term1 + term2 + term3

        return thrust_to_weight

    @staticmethod
    def constant_velocity_turn_imp(wing_loading: float | ndarray, c_d_min: float, k: float, n: float, v_cruise: float, rho: float) -> float | ndarray:
        q = Physics.calc_dynamic_pressure(rho, v_cruise)

        thrust_to_weight = q * ((c_d_min / wing_loading) + k * pow(n / q, 2) * wing_loading)

        return thrust_to_weight

    @staticmethod
    def cruise_imp(wing_loading: float | ndarray, c_d_min: float, k: float, v_cruise: float, rho: float) -> float | ndarray:
        q = Physics.calc_dynamic_pressure(rho, v_cruise)

        thrust_to_weight = q * c_d_min * (1 / wing_loading) + k * (1 / q) * wing_loading

        return thrust_to_weight

    @staticmethod
    def service_ceiling_imp(wing_loading: float | ndarray, c_d_min: float, k: float, v_upsilon_kcas: float, altitude_m: float) -> float | ndarray:
        return AircraftPerformance.rate_of_climb_imp(wing_loading, 1.667, c_d_min, k, v_upsilon_kcas, altitude_m)

    @staticmethod
    def landing_distance_imp(target_landing_distance_ft: float, airport_altitude_m: float, cl_max, cl_landing, cd_landing: float, mu: float = 0.03, h_obst_ft: float = 35, free_roll_time_s: float = 1, thrust_loading_ground: float = 0.0) -> float:
        rho = ISAtmosphere().get_density(airport_altitude_m)
        rho_imp = Convert().convert(Convert().convert(rho, Convert.KG, Convert.SLUG), Convert.METER, Convert.FOOT, -3)

        a = rho_imp * cl_max

        wing_loading = 0.0
        actual_landing_distance_ft = 1

        while actual_landing_distance_ft <= target_landing_distance_ft:
            wing_loading = wing_loading + 0.01
            term1 = 19.08 * h_obst_ft
            term2 = 0.007923 + 1.556 * free_roll_time_s * np.sqrt(a / wing_loading)
            term3 = 1.21 / (PhysicalConstant.G * ((0.605 / cl_max) * (cd_landing - mu * cl_landing) + mu - thrust_loading_ground))

            actual_landing_distance_ft = term1 + (term2 + term3) * (wing_loading / a)

        wing_loading = wing_loading - 0.01

        return wing_loading
