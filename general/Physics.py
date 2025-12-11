from numpy import ndarray


class Physics:
    @staticmethod
    def calc_dynamic_pressure(rho: float | ndarray, v_ft_per_s: float | ndarray) -> float | ndarray:
        q = 0.5 * rho * pow(v_ft_per_s, 2)
        return q