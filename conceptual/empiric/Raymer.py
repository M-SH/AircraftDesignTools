class Raymer:
    @staticmethod
    def calc_oswald_factor(ar: float) -> float:
        return 1.78 * (1 - 0.045 * pow(ar, 0.68)) - 0.64
