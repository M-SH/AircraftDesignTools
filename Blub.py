class System:
    def __init__(self, name, position=(0, 0, 0)):
        """
        Initialisiert ein System mit einem Namen und einer absoluten Position im 3D-Raum.

        :param name: Name des Systems
        :param position: Absolute Position des Systems (x, y, z)
        """
        self.name = name
        self.position = position
        self.subsystems = []

    def add_subsystem(self, subsystem, relative_position=(0, 0, 0)):
        """
        Fügt ein Sub-System mit einer relativen Position im 3D-Raum hinzu.

        :param subsystem: Instanz des Sub-Systems
        :param relative_position: Relative Position des Sub-Systems (x, y, z)
        """
        subsystem.relative_position = relative_position
        self.subsystems.append(subsystem)

    def absolute_position(self):
        """
        Berechnet die absolute Position jedes Sub-Systems, einschließlich
        aller hierarchisch untergeordneten Systeme im 3D-Raum.

        :return: Dictionary mit den absoluten Positionen der Sub-Systeme
        """
        absolute_positions = {}
        for subsystem in self.subsystems:
            x_abs = self.position[0] + subsystem.relative_position[0]
            y_abs = self.position[1] + subsystem.relative_position[1]
            z_abs = self.position[2] + subsystem.relative_position[2]
            absolute_positions[subsystem.name] = (x_abs, y_abs, z_abs)
            # Berechne auch die Positionen der Subsysteme dieses Subsystems
            if subsystem.subsystems:
                nested_positions = subsystem.absolute_position()
                for name, pos in nested_positions.items():
                    absolute_positions[f"{subsystem.name} -> {name}"] = pos
        return absolute_positions


# Beispiel der Verwendung
if __name__ == "__main__":
    main_system = System("Hauptsystem", (10, 20, 30))

    subsystem1 = System("Sub-System 1", (0, 0, 0))
    subsystem2 = System("Sub-System 2", (0, 0, 0))

    # Füge Subsysteme zu Sub-System 1 hinzu
    nested_subsystem = System("Unter-Subsystem 1", (0, 0, 0))
    subsystem1.add_subsystem(nested_subsystem, (5, 5, 5))

    # Füge Subsysteme zum Hauptsystem hinzu
    main_system.add_subsystem(subsystem1, (5, 5, 5))
    main_system.add_subsystem(subsystem2, (-3, 2, 1))

    abs_positions = main_system.absolute_position()
    print(abs_positions)