from numpy import ndarray
from plotly.graph_objs import Figure
from plotly.graph_objs.layout import Shape

class CabinArtifact:
    def __init__(self, name: str = "Artifact", position: ndarray[float] = None, label: str = None, width: float = 1.0, length: float = 1.0):
        if position is None:
            position = [0, 0]

        if label is None:
            label = name

        self.name = name
        self.label = label
        self.position = position
        self.width = width
        self.length = length

    def draw(self, fig: Figure):
        s = Shape(type="rect", x0=self.position[0], y0=self.position[1], x1=self.position[0] + self.length, y1 = self.position[1] + self.width,
                      line=dict(color="black", width = 2), label=dict(text=self.label))
        fig.add_shape(s)
