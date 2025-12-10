import numpy as np
from matplotlib import pyplot as plt
from numpy import ndarray, linspace

from general.Convert import Convert
from general.ISAtmosphere import ISAtmosphere
from general.Singleton import Singleton

class Airspeeds(Singleton):
    """
    Contains static methods to convert between air speeds (KCAS, KEAS, KTAS).
    """
    @staticmethod
    def to_keas(v_kcas: float | ndarray, altitude_m: float | ndarray) -> float | ndarray:
        """Converts knots **calibrated** air speed to knots **equivalent** air speed for a given altitude.

        *Reference: Gudmundsson, General Aviation Aircraft Design, 2nd Edition, Equation (17-23), Page 763*

        :param v_kcas: Flight speed in knots calibrated airspeed (KCAS).
        :param altitude_m: Flight altitude in meters.
        :return: Flight speed in knots equivalent airspeed (KEAS).
        """
        q_c = Airspeeds.calc_compr_dynamic_pressure(v_kcas, altitude_m)
        p = ISAtmosphere().get_pressure(altitude_m)

        term1 = v_kcas * np.sqrt(ISAtmosphere().get_pressure_ratio(altitude_m))
        term2 = pow((q_c / p) + 1, 0.2857) - 1
        term3 = pow((q_c / ISAtmosphere().p0) + 1, 0.2857) - 1

        v_keas = term1 * np.sqrt(term2 / term3)

        return v_keas

    @staticmethod
    def to_ktas(v_keas: float | ndarray, altitude_m: float | ndarray) -> float | ndarray:
        """Converts knots **equivalent** air speed to knots **true** air speed for a given altitude.

        *Reference: Gudmundsson, General Aviation Aircraft Design, 2nd Edition, Equation (17-30), Page 764*

        :param v_keas: Flight speed in knots equivalent airspeed (KEAS).
        :param altitude_m: Flight altitude in meters.
        :return: Flight speed in knots true airspeed (KTAS).
        """
        return v_keas / np.sqrt(ISAtmosphere().get_density_ratio(altitude_m))

    @staticmethod
    def calc_compr_dynamic_pressure(v_kcas: float, altitude_m: float | ndarray = 0.0) -> float | ndarray:
        """ Calculates compressible dynamic pressure for a given calibrated flight speed and altitude.

        *Reference: Gudmundsson, General Aviation Aircraft Design, 2nd Edition, Equation (17-25), Page 763*

        :param v_kcas: Flight speed in knots true airspeed (KTAS).
        :param altitude_m: Flight altitude in meters.
        :return: Compressible dynamic pressure q_c in ???
        """
        p = ISAtmosphere().get_pressure(altitude_m)
        ma = Airspeeds.calc_mach_from_kcas(v_kcas, altitude_m)

        q_c = p * (pow(1 + 0.2 * pow(ma, 2), 3.5) - 1)

        return q_c

    @staticmethod
    def calc_mach_from_ktas(v_ktas: float | ndarray, altitude_m: float | ndarray = 0.0) -> float | ndarray:
        """ Calculates the Mach number for a given **true** airspeed speed and altitude.

        *Reference: Gudmundsson, General Aviation Aircraft Design, 2nd Edition, Equation (17-26), Page 763*

        :param v_ktas: Flight speed in knots true airspeed (KTAS).
        :param altitude_m: Flight altitude in meters.
        :return: Mach number.
        """
        return v_ktas / ISAtmosphere().get_speed_of_sound(altitude_m)

    @staticmethod
    def calc_mach_from_kcas(v_kcas: float | ndarray, altitude_m: float | ndarray = 0.0) -> float | ndarray:
        """ Calculates the Mach number for a given **calibrated** airspeed speed and altitude.

        *Reference: Gudmundsson, General Aviation Aircraft Design, 2nd Edition, Equation (17-28), Page 764*

        :param v_kcas: Flight speed in knots calibrated airspeed (KTAS).
        :param altitude_m: Flight altitude in meters.
        :return: Mach number.
        """
        pressure_ratio = ISAtmosphere().get_pressure_ratio(altitude_m)
        return 2.236 * np.sqrt(pow(((pow(1 + 4.575 * pow(10, -7) * pow(v_kcas, 2), 3.5) - 1) / pressure_ratio) + 1, 0.2857) - 1)

    @staticmethod
    def plot_vs_altitude(v_kcas: float = 265.0) -> None:
        """Creates a plot showing KCAS, KEAS and KTAS vs. altitude.

        If no speed is specified, 265 KCAS is used to create the plot. (This reproduces Figure 17.10 in Gudmundsson.)

        *Reference: Gudmundsson, General Aviation Aircraft Design, 2nd Edition, Figure (17.10), Page 765*

        :param v_kcas: Flight speed in knots calibrated airspeed (KTAS).
        """
        altitudes = linspace(0,50000, 100) # ft
        altitudes_m = Convert().convert(altitudes, Convert.FOOT, Convert.METER)

        v_kcas = np.array([v_kcas] * 100)
        v_keas = Airspeeds.to_keas(v_kcas, altitudes_m)
        v_ktas = Airspeeds.to_ktas(v_keas, altitudes_m)

        plt.figure()
        plt.plot(v_kcas, altitudes, label='CAS')
        plt.plot(v_keas, altitudes, label='EAS')
        plt.plot(v_ktas, altitudes, label='TAS')
        plt.xlabel('Airspeed [knots]')
        plt.ylabel('Altitude [1000ft]')
        plt.title('Airspeeds vs. Altitude')
        plt.grid(True)
        plt.xlim(0.0, 700.0)
        plt.ylim(0.0, 50000.0)
        plt.legend(loc='upper right')
        plt.show()

if __name__ == "__main__":
    Airspeeds.plot_vs_altitude()
