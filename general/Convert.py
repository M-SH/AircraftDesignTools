
from numpy import ndarray

from general.Singleton import Singleton
from general.StringConstant import StringConstant


class Convert(metaclass=Singleton):
    METER = "m"
    FOOT = "ft"
    KG = "kg"
    POUND = "lbf"
    SLUG = "slug"

    METER_TO_FEET = 3.281
    KG_TO_POUND = 2.20462
    KG_TO_SLUG = 0.0685218


    def __init__(self):
        self.conversion_dict = dict[str, float]()
        self.conversion_dict[Convert._create_key(Convert.METER, Convert.FOOT)] = Convert.METER_TO_FEET
        self.conversion_dict[Convert._create_key(Convert.FOOT, Convert.METER)] = 1.0 / Convert.METER_TO_FEET
        self.conversion_dict[Convert._create_key(Convert.KG, Convert.POUND)] = Convert.KG_TO_POUND
        self.conversion_dict[Convert._create_key(Convert.POUND, Convert.KG)] = 1.0 / Convert.KG_TO_POUND
        self.conversion_dict[Convert._create_key(Convert.KG, Convert.SLUG)] = Convert.KG_TO_SLUG
        self.conversion_dict[Convert._create_key(Convert.SLUG, Convert.KG)] = 1.0 / Convert.KG_TO_SLUG

    def convert(self, value: float | ndarray, current_unit: str, target_unit: str, dim: int = 1) -> float | ndarray:
        conversion_factor = self.conversion_dict.get(self._create_key(current_unit, target_unit))
        return value * pow(conversion_factor, float(dim))

    @staticmethod
    def _create_key(current_unit: str, target_unit: str):
        return current_unit + StringConstant.UNDERSCORE + target_unit