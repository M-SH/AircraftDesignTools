from abc import ABC

import numpy as np

from model.Assembly import Assembly


class System(Assembly):
    def __init__(self, name: str):
        super().__init__(name)
        self.mass_kg = 0.0

    def get_mass_kg(self) -> float:
        return self.mass_kg
