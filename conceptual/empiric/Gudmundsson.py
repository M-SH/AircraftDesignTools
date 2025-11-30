from numpy import ndarray


class Gudmundsson:
    @staticmethod
    def calc_best_roc_speed_single_engine_cas(w_s: float | ndarray) -> float | ndarray:
        # gudmundsson 2nd edition equ. 3.3
        v_upsilon = 43.591 + 2.2452 * w_s
        return v_upsilon

    @staticmethod
    def calc_best_roc_speed_twin_engine_cas(w_s: float | ndarray) -> float | ndarray:
        # gudmundsson 2nd edition equ. 3.4
        v_upsilon = 69.952 + 1.3402 * w_s
        return v_upsilon

    @staticmethod
    def calc_best_roc_speed_business_jet_cas(w_s: float | ndarray) -> float | ndarray:
        # gudmundsson 2nd edition equ. 3.5
        v_upsilon = 79.016 + 1.2722 * w_s
        return v_upsilon