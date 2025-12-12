import numpy as np
import plotly.graph_objects as go
from numpy import ndarray
from plotly.graph_objs import Figure

from general.Convert import Convert
from model.Assembly import Assembly
from model.aircraft.CabinArtifact import CabinArtifact
from model.aircraft.Door import Door
from model.aircraft.Seat import Seat


class Cabin(Assembly):
    def __init__(self, name: str = "Cabin", width: float = 2.5, length: float = 16.5):
        super().__init__(name)
        self.width = width
        self.length = length
        self.doors = []
        self.artifacts = []
        self.seats = []

    def generate_seats(self, position, seat_pitch_inch: float = 29, abreast: int = 2, rows: int = 17):
        seat_pitch_m = Convert().convert(seat_pitch_inch, Convert.INCH, Convert.METER)
        current_position = position
        for i in range(0, rows):
            for j in range(0, abreast):
                seat = Seat(current_position)
                self.seats.append(seat)
                current_position = [current_position[0], current_position[1] + seat.width]

            current_position = [current_position[0] + seat_pitch_m, position[1]]

    def add_door(self, name: str = "Door", side: str = "left", x: float = 1.0, width: float = 1.0, height: float = 1.0):
        y = 0.0
        if side == "right":
            y = self.width

        self.doors.append(Door(name, position=[x, y], width=width, height=height))

    def draw(self, fig: Figure = None):
        if fig is None:
            fig = go.Figure()

        fig.add_shape(type="rect", x0=0, y0=0, x1=self.length, y1=self.width, line=dict(color="black", width=2))

        for seat in self.seats:
            seat.draw(fig)

        for artifact in self.artifacts:
            artifact.draw(fig)

        for door in self.doors:
            door.draw(fig)

        axis_limit_x = np.ceil(self.length)
        axis_ticks_x = int((axis_limit_x + 1) * 4)

        fig.update_xaxes(
            range=[-1, axis_limit_x],
            nticks=axis_ticks_x
        )

        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
            range=[-1, 3],
            nticks = 20
        )

        fig.update_layout(newshape=dict(line_color='red', opacity=1, label=dict(texttemplate="%{length:.3f}m")))
        fig.update_layout(hovermode="x")

        fig.write_html('first_figure.html', config={'modeBarButtonsToAdd':['drawline', 'eraseshape']}, auto_open=True)


if __name__ == "__main__":
    cabin_atr72_68seats = Cabin("ATR72-600 68PAX", width=2.57, length=16.52)
    cabin_atr72_68seats.artifacts.append(CabinArtifact("Front Cargo Compartment Left", [0.0, 0.0], width=0.892, length=1.916))
    cabin_atr72_68seats.artifacts.append(CabinArtifact("Front Cargo Compartment Right", [0.0, 1.678], width=0.892, length=1.916))
    cabin_atr72_68seats.artifacts.append(
        CabinArtifact("Galley Left",[15.339, 0.0], label="G", width=0.988, length=0.301))
    cabin_atr72_68seats.artifacts.append(
        CabinArtifact("Galley Right",[15.339, 2.57 - 0.988], label="G", width=0.988, length=0.482))
    cabin_atr72_68seats.generate_seats([2.55, 0.1], seat_pitch_inch=31, abreast=2, rows=17)
    cabin_atr72_68seats.generate_seats([2.55, 1.5], seat_pitch_inch=31, abreast=2, rows=17)
    cabin_atr72_68seats.add_door("Type III Exit", "left", x=2.313, width=0.51, height=0.91)
    cabin_atr72_68seats.add_door("Type III Exit", "right", x=2.313, width=0.51, height=0.91)
    cabin_atr72_68seats.add_door("Cargo Door", "left", x=0.361, width=1.275, height=1.53)
    cabin_atr72_68seats.add_door("Entrance Door", "left", x=15.724, width=0.75, height=1.75)
    cabin_atr72_68seats.add_door("Service Door", "right", x=15.869, width=0.61, height=1.22)

    cabin_atr72_68seats.draw()