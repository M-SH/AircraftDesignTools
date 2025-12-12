from numpy import ndarray
from plotly.graph_objs import Figure
from plotly.graph_objs.layout import Shape

from general.StringConstant import StringConstant


class Door:
    def __init__(self, name: str = "Door", position: ndarray[float] = None, width: float = 1.0, height: float = 1.0):
        if position is None:
            position = [0, 0]

        self.name = name
        self.label = name
        self.position = position
        self.width = width
        self.height = height

    def draw(self, fig: Figure):
        text_pos = "bottom center"
        if self.position[1] == 0.0:
            text_pos = "top center"

        svg_path = "M " + str(self.position[0]) + StringConstant.BLANK + str(self.position[1]) + " L " + str(self.position[0] + self.width) + StringConstant.BLANK + str(self.position[1])

        s = Shape(type="path", path=svg_path,
                  line=dict(color="white", width=1.5), label=dict(text=self.label, padding=20, textposition=text_pos))
        fig.add_shape(s)