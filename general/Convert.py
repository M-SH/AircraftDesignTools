
from numpy import ndarray

from general.Singleton import Singleton
from general.StringConstant import StringConstant


class Convert(metaclass=Singleton):
    METER = "m"
    FOOT = "ft"
    INCH = "inch"

    KG = "kg"
    POUND = "lbf"
    SLUG = "slug"

    KNOT = "kts"
    FT_PER_S = "ft_per_s"
    M_PER_S = "m_per_s"

    METER_TO_FEET = 3.281
    METER_TO_INCH = 39.37
    KG_TO_POUND = 2.20462
    KG_TO_SLUG = 0.0685218
    KTS_TO_FT_PER_S = 1.68781


    def __init__(self):
        self.conversion_dict = dict[str, float]()
        self.conversion_dict[Convert._create_key(Convert.METER, Convert.FOOT)] = Convert.METER_TO_FEET
        self.conversion_dict[Convert._create_key(Convert.FOOT, Convert.METER)] = 1.0 / Convert.METER_TO_FEET
        self.conversion_dict[Convert._create_key(Convert.METER, Convert.INCH)] = Convert.METER_TO_INCH
        self.conversion_dict[Convert._create_key(Convert.INCH, Convert.METER)] = 1.0 / Convert.METER_TO_INCH
        self.conversion_dict[Convert._create_key(Convert.KG, Convert.POUND)] = Convert.KG_TO_POUND
        self.conversion_dict[Convert._create_key(Convert.POUND, Convert.KG)] = 1.0 / Convert.KG_TO_POUND
        self.conversion_dict[Convert._create_key(Convert.KG, Convert.SLUG)] = Convert.KG_TO_SLUG
        self.conversion_dict[Convert._create_key(Convert.SLUG, Convert.KG)] = 1.0 / Convert.KG_TO_SLUG
        self.conversion_dict[Convert._create_key(Convert.KNOT, Convert.FT_PER_S)] = Convert.KTS_TO_FT_PER_S
        self.conversion_dict[Convert._create_key(Convert.FT_PER_S, Convert.KNOT)] = 1.0 / Convert.KTS_TO_FT_PER_S
        self.conversion_dict[Convert._create_key(Convert.M_PER_S, Convert.FT_PER_S)] = Convert.METER_TO_FEET
        self.conversion_dict[Convert._create_key(Convert.FT_PER_S, Convert.M_PER_S)] = 1.0 / Convert.METER_TO_FEET

    def convert(self, value: float | ndarray, current_unit: str, target_unit: str, dim: int = 1) -> float | ndarray:
        conversion_factor = self.conversion_dict.get(self._create_key(current_unit, target_unit))
        return value * pow(conversion_factor, float(dim))

    @staticmethod
    def _create_key(current_unit: str, target_unit: str):
        return current_unit + StringConstant.UNDERSCORE + target_unit