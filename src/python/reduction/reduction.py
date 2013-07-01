
from numpy import dot, diag, floor, trace, eye, linspace
from numpy.linalg import svd, pinv, inv, norm

from models.signal_model import ProbabilisticSignalModel, FuzzySignalModel

class Reduction(object):
    def __init__(self, scheme, *args, **kwargs):
        self.scheme = scheme

    def compute(self):
        raise NotImplementedError()

class LinearReduction(Reduction):
    def __init__(self, scheme, ideal_operator, *args, **kwargs):
        super(LinearReduction, self).__init__(scheme, *args, **kwargs)
        self.ideal_operator = ideal_operator
        self.reduction_operator = None

    def compute(self):
        raise NotImplementedError()

class ProbabilisticLinearReduction(LinearReduction):
    def __init__(self, scheme, ideal_operator, *args, **kwargs):
        super(ProbabilisticLinearReduction, self).__init__(
            scheme, ideal_operator, *args, **kwargs)

    def _noise_covariance(self):
        return self.scheme.noise.cov_matrix

    def prepare(self):
        """Numerically stable way to compute reduction.
        First compute SVD to look at the distribution
        of singular values. Then null the smallest singular values.
        """
        operator = self.scheme.operator
        noise_covariance = self._noise_covariance()

        supposed_invertable = dot(dot(operator.T, inv(noise_covariance)), operator)

        tolerance = 1e-8
        right_singular, singular_values, left_singular = svd(supposed_invertable)
        inverse_singular_values = diag(pinv(diag(singular_values), tolerance))
        max_inverse_value_index = floor(2. * (len(inverse_singular_values) / 3.))
        truncated_singular_values_indices = \
            inverse_singular_values > inverse_singular_values[max_inverse_value_index]
        inverse_singular_values[truncated_singular_values_indices] = 0
        truncated_operator_inverse = dot(dot(right_singular, diag(inverse_singular_values)), left_singular)

        self.reduction_operator = dot(dot(self.ideal_operator, truncated_operator_inverse), operator.T)

    def compute(self):
        if None == self.reduction_operator:
            self.prepare()

        measurement = self.scheme.measurement.value
        self.estimate = ProbabilisticSignalModel(value=dot(self.reduction_operator, measurement))

        signal = self.scheme.signal.value
        dimensions, = (signal.shape)

        ideal_estimate = dot(self.ideal_operator, signal)
        error_value = dot(ideal_estimate - self.estimate.value, ideal_estimate - self.estimate.value)
        self.error = ProbabilisticSignalModel(value=error_value / dimensions)

class FuzzyLinearReduction(ProbabilisticLinearReduction):
    def __init__(self, scheme, ideal_operator, *args, **kwargs):
        super(FuzzyLinearReduction, self).__init__(
            scheme, ideal_operator, *args, **kwargs)

    def _noise_covariance(self):
        return self.scheme.noise.coordinated_model.cov_matrix

    def compute(self):
        if None == self.reduction_operator:
            self.prepare()

        measurement = self.scheme.measurement.value
        self.estimate = FuzzySignalModel(value=dot(self.reduction_operator, measurement))

        signal = self.scheme.signal.value
        dimensions, = (signal.shape)

        ideal_estimate = dot(self.ideal_operator, signal)
        error_value = dot(ideal_estimate - self.estimate.value, ideal_estimate - self.estimate.value)
        self.error = FuzzySignalModel(value=error_value / dimensions)

class LinearRegularization(LinearReduction):
    def __init__(self, scheme, ideal_operator, *args, **kwargs):
        super(LinearRegularization, self).__init__(
            scheme, ideal_operator, *args, **kwargs)
        self.regularizer = 0.

    def do_compute_with_regularizer(self, regularizer):
        operator = self.scheme.operator
        measurement = self.scheme.measurement.value

        supposed_invertable = dot(operator.T, operator)
        supposed_invertable += regularizer * eye(supposed_invertable.shape[0])

        tolerance = 1e-8
        right_singular, singular_values, left_singular = svd(supposed_invertable)
        inverse_singular_values = diag(pinv(diag(singular_values), tolerance))
        max_inverse_value_index = floor(2. * (len(inverse_singular_values) / 3.))
        truncated_singular_values_indices = \
            inverse_singular_values > inverse_singular_values[max_inverse_value_index]
        inverse_singular_values[truncated_singular_values_indices] = 0
        truncated_operator_inverse = dot(dot(right_singular, diag(inverse_singular_values)), left_singular)

        self.reduction_operator = dot(dot(self.ideal_operator, truncated_operator_inverse), operator.T)

        self.estimate = ProbabilisticSignalModel(value=dot(self.reduction_operator, measurement))

        signal = self.scheme.signal.value
        dimensions, = (signal.shape)

        ideal_estimate = dot(self.ideal_operator, signal)
        error_value = dot(ideal_estimate - self.estimate.value, ideal_estimate - self.estimate.value)
        self.error = ProbabilisticSignalModel(value=error_value / dimensions)
        self.regularizer = regularizer

    def compute(self):
        left = 0.
        right = 10.
        for iteration in xrange(10):
            center = left + (right - left) / 2.

            self.do_compute_with_regularizer(left)
            left_error = self.error.value

            self.do_compute_with_regularizer(center)
            center_error = self.error.value

            self.do_compute_with_regularizer(right)
            right_error = self.error.value

            if left_error >= center_error >= right_error:
                left = center
            elif left_error <= center_error <= right_error:
                right = center
            else:
                left_center = left + (center - left) / 2.
                self.do_compute_with_regularizer(left_center)
                left_center_error = self.error.value

                if left_error > left_center_error > center_error:
                    left = center
                else:
                    right = center
