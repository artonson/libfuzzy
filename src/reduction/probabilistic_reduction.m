function [estimate, error_estimate] = probabilistic_reduction(measurement, ...
        operator, varargin)
%PROBABILISTIC_REDUCTION Compute linear estimate via reduction
%   PROBABILISTIC_REDUCTION computes the best linear estimate for measurement
%   given the measurement operator, (optional) measurement expected value 
%   and covariance operator and (optional) noise covariance matrix. 

numvarargs = length(varargin);
if numvarargs > 3
    error('libfuzzy:probabilistic_reduction:too_many_inputs', ...
        'requires at most 3 optional inputs')
    
elseif numvarargs == 1
    % Compute a more accurate estimate using the noise covariance matrix.
    noise_covariance = varargin{1};
    estimate = pinv(operator' * inv(noise_covariance) * operator) *...
        operator' * measurement;
    error_estimate = trace(operator' * inv(noise_covariance) * operator);
    
elseif numvarargs == 2
    % Compute a more accurate estimate using the noise covariance matrix
    % and the expected value.
    noise_covariance = varargin{1};
    signal_covariance = varargin{2};
    reduction_operator = signal_covariance * operator' * ...
        inv(operator * signal_covariance * operator' + noise_covariance);
    estimate = reduction_operator * measurement;
    error_estimate = trace(signal_covariance - ...
        reduction_operator * operator * signal_covariance);
    
elseif numvarargs == 3
    % Compute a more accurate estimate using the noise covariance matrix,
    % the expected value, and the signal covariance matrix.
    noise_covariance = varargin{1};
    signal_covariance = varargin{2};
    expected_signal = varargin{3};
    reduction_operator = signal_covariance * operator' * ...
        inv(operator * signal_covariance * operator' + noise_covariance);
    estimate = expected_signal + reduction_operator * ...
        (measurement - operator * expected_signal);
    error_estimate = trace(signal_covariance - ...
        reduction_operator * operator * signal_covariance);

else 
    % Compute basic least squares estimate.
    estimate = pinv(operator) * measurement;
    error_estimate = trace(operator' * operator);
    
end

end
