from numpy import ndarray

from general.Convert import Convert
from general.PhysicalConstant import PhysicalConstant

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

        term1 = (1.21 / (g_imp * rho * cl_max * ground_run_distance)) ** wing_loading
        term2 = (0.605 / cl_max) * (cd_to - my * cl_to)
        term3 = my

        thrust_to_weight = term1 + term2 + term3

        return thrust_to_weight
