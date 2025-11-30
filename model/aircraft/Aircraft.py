from model.Assembly import Assembly


class Aircraft(Assembly):
    def __init__(self, name: str):
        super().__init__(name)
        self.ar = 5
        self.s_ref = 120
