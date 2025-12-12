from model.aircraft.CabinArtifact import CabinArtifact

class Seat(CabinArtifact):
    def __init__(self, position=None, width: float = 0.46, length: float = 0.46, label: str = None):
        super().__init__("Seat", position, width=width, length=length, label=label)
        # self.label = "S"