function [estimate, error_estimate] = probabilistic_reduction(measurement, ...
        operator, ideal_operator, signal_model, noise_model)
%PROBABILISTIC_REDUCTION Compute appropriate linear estimate via reduction
%   PROBABILISTIC_REDUCTION computes the best linear estimate for measurement
%   given the measurement operator, (optional) measurement expected value 
%   and covariance operator and (optional) noise covariance matrix. 


args = {};

if isfield(noise_model, 'covariance')
    args{1} = noise_model.covariance;
    if isfield(signal_model, 'covariance')
        args{2} = signal_model.covariance;
        if isfield(signal_model, 'expectation')
            args{3} = signal_model.expectation;
        end
    end
end

[estimate, error_estimate] = do_probabilistic_reduction(measurement, ...
        operator, ideal_operator, args);

    function [estimate, error_estimate] = do_probabilistic_reduction(measurement, ...
        operator, ideal_operator, args)
%DO_PROBABILISTIC_REDUCTION Compute appropriate linear estimate via reduction
%   DO_PROBABILISTIC_REDUCTION computes the best linear estimate for measurement
%   given the measurement operator, (optional) measurement expected value 
%   and covariance operator and (optional) noise covariance matrix. 

    numvarargs = length(args);
    if numvarargs > 3
        error('libfuzzy:probabilistic_reduction:too_many_inputs', ...
            'requires at most 3 optional inputs')

    elseif numvarargs == 1
        % Compute a more accurate estimate using the noise covariance matrix.
        noise_covariance = args{1};
        reduction_operator = ideal_operator * ...
            pinv(operator' * inv(noise_covariance) * operator, 1e-5) *...
            operator';
        estimate = reduction_operator * measurement;
        error_estimate = trace(operator' * inv(noise_covariance) * operator);

    elseif numvarargs == 2
        % Compute a more accurate estimate using the noise covariance matrix
        % and the expected value.
        noise_covariance = args{1};
        signal_covariance = args{2};
        reduction_operator = signal_covariance * operator' * ...
            pinv(operator * signal_covariance * operator' + noise_covariance, 1e-8);
        estimate = reduction_operator * measurement;
        error_estimate = trace(signal_covariance - ...
            reduction_operator * operator * signal_covariance);

    elseif numvarargs == 3
        % Compute a more accurate estimate using the noise covariance matrix,
        % the expected value, and the signal covariance matrix.
        noise_covariance = args{1};
        signal_covariance = args{2};
        expected_signal = args{3};
        reduction_operator = signal_covariance * operator' * ...
            pinv(operator * signal_covariance * operator' + noise_covariance, 1e-8);
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

end
