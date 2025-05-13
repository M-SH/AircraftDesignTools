from model.System import System
from model.aircraft.Aircraft import Aircraft

if __name__ == "__main__":
    do328 = Aircraft("D328eco")
    do328.set_position_m(0, 0, 0)
    airframe = System("Airframe")
    airframe.mass_kg = 5000.0
    airframe.set_position_m(1, 0, 0)
    airframe.set_rotation_deg(45, 0, 0)
    engines = System("Engines")
    engines.mass_kg = 1000.0
    engines.set_position_m(0, 2, 0)
    payload = System("Payload")
    payload.mass_kg = 4200.0
    fuel = System("Fuel")
    fuel.mass_kg = 900.0

    do328.add_sub_system(airframe)
    airframe.add_sub_system(engines)
    do328.add_sub_system(payload)
    do328.add_sub_system(fuel)

    print(do328.get_mass_kg())
    print(engines.get_absolute_position_m())
