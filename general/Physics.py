from numpy import ndarray


class Physics:
    @staticmethod
    def calc_dynamic_pressure(rho: float | ndarray, v: float | ndarray) -> float | ndarray:
        q = 0.5 * rho * pow(v, 2)
        return q