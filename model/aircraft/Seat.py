from plotly.graph_objs import Figure
import plotly.graph_objects as go

class Seat:
    def __init__(self, position = None, width: float = 0.4572, length: float = 0.4572):
        if position is None:
            position = [0, 0]
        self.position = position
        self.width = width
        self.length = length

    def draw(self, fig: Figure):
        fig.add_trace(go.Scatter(x=[self.position[0], self.position[0], self.position[0] + self.length, self.position[0] + self.length], y=[self.position[1], self.position[1] + self.width, self.position[1] + self.width, self.position[1]], fill="toself"))
