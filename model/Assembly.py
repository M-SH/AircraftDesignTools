from abc import ABC

import numpy as np


class Assembly(ABC):
    def __init__(self, name: str):
        self.name = name
        self.position_relative_m = np.array([0.0, 0.0, 0.0])
        self.rotation_relative_m = np.array([0.0, 0.0, 0.0])

        self.parent_system = None
        self.sub_systems = []

    def set_position_m(self, x: float, y: float, z: float = 0.0):
        self.position_relative_m = np.array([x, y, z])

    def set_rotation_deg(self, rot_x: float, rot_y: float, rot_z: float = 0.0):
        self.position_relative_m = np.array([rot_x, rot_y, rot_z])

    def get_absolute_position_m(self) -> np.array:
        if self.parent_system is not None:
            return self.parent_system.get_absolute_position_m() + self.position_relative_m
        else:
            return self.position_relative_m

    def get_mass_kg(self) -> float:
        return sum([s.mass_kg for s in self.sub_systems])

    def add_sub_system(self, sub_system: 'Assembly'):
        self.sub_systems.append(sub_system)
        sub_system.parent_system = self

