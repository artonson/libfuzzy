import numpy as np

class MeasurementScheme:
    def measure(self):
        raise NotImplementedError()

class LinearMeasurementScheme(MeasurementScheme):
    def __init__(self, operator, signal, noise):
        self.operator = operator
        self.signal = signal
        self.noise = noise
        self.measurement = None

    def measure(self):
        signal_class = self.signal.__class__
        self.signal.sample()
        self.noise.sample()
        measurement_value = np.dot(self.operator, self.signal.value) + self.noise.value
        self.measurement = signal_class(value=measurement_value)
        return self.measurement, self.signal

