import plotly.graph_objects as go

from model.aircraft.Seat import Seat


class Cabin:
    def __init__(self, width: float = 2.5, length: float = 16.5):
        self.width = width
        self.length = length
        self.seats = []

    def generate_seats(self, position, seat_pitch_inch: float = 29, abreast: int = 2, rows: int = 18):
        seat_pitch_m = seat_pitch_inch * 0.0254
        current_position = position
        for i in range(0, rows):
            for j in range(0, abreast):
                seat = Seat(current_position)
                self.seats.append(seat)
                current_position = [current_position[0], current_position[1] + seat.width]

            current_position = [current_position[0] + seat_pitch_m, position[1]]

    def draw(self):
        fig = go.Figure(go.Scatter(x=[0, 0, self.length, self.length], y=[0, self.width, self.width, 0], fill="toself"))

        for seat in self.seats:
            seat.draw(fig)

        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
        )

        fig.show()

if __name__ == "__main__":
    cabin = Cabin()
    cabin.generate_seats([1.0, 0.1])
    cabin.generate_seats([1.0, 1.5])
    cabin.draw()