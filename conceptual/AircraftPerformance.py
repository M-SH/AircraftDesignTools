from numpy import ndarray

from general.Convert import Convert
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
    def rate_of_climb_imp(wing_loading: float | ndarray, roc: float, c_d_min: float, k: float, v_upsilon: float, rho: float) -> float | ndarray:
        q = Physics.calc_dynamic_pressure(rho, v_upsilon)

        term1 = roc / v_upsilon
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
    def service_ceiling_imp(wing_loading: float | ndarray, c_d_min: float, k: float, v_upsilon: float, rho: float) -> float | ndarray:
        return AircraftPerformance.rate_of_climb_imp(wing_loading, 1.667, c_d_min, k, v_upsilon, rho)
