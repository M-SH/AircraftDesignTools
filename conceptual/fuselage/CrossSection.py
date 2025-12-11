import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse


class CrossSection:
    def __init__(self):
        self.cabin_floor_width_m: float = 0.0

        self.bubble_upper = CrossSection.Bubble([0.0, 1.0], 1.2, 1.0)
        self.bubble_lower = CrossSection.Bubble([0.0, 0.0], 1.0, 0.5)

        self.width_m: float = 2.0
        self.height_upper_m: float = 2.0
        self.height_lower_m: float = 2.0

    def draw(self):
        [x_upper, y_upper] = self.bubble_upper.calculate_plot_coordinates()
        [x_lower, y_lower] = self.bubble_lower.calculate_plot_coordinates()

        plt.plot(x_upper, y_upper)
        plt.plot(x_lower, y_lower)

        #u = 1.  # x-position of the center
        #v = 0.5  # y-position of the center
        #a = 2.  # radius on the x-axis
        #b = 1.5  # radius on the y-axis

        #t = np.linspace(0, 2 * np.pi, 100)
        #plt.plot(u + a * np.cos(t), v + b * np.sin(t))
        plt.grid(color='lightgray', linestyle='--')
        plt.show()

    class Bubble:
        def __init__(self, center = None, width: float = 1.0, height: float = 1.0):
            if center is None:
                center = [0.0, 0.0]

            self.center = center
            self.width = width
            self.height = height

        def calculate_plot_coordinates(self):
            t = np.linspace(0, 2 * np.pi, 100)
            x = self.center[0] + (self.width / 2.0) * np.cos(t)
            y = self.center[1] + (self.height / 2.0) * np.sin(t)
            return [x, y]
