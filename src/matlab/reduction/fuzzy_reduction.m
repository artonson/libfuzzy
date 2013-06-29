function [estimate, error_estimate] = fuzzy_reduction(measurement, ...
        operator, ideal_operator, signal_model, noise_model)
%FUZZY_REDUCTION Summary of this function goes here
%   Detailed explanation goes here

args = {};

probabilistic_noise_model = noise_model.coordinated_probabilistic_model;
probabilistic_signal_model = signal_model.coordinated_probabilistic_model;

if isfield(probabilistic_noise_model, 'covariance')
    args{1} = probabilistic_noise_model.covariance;
    if isfield(probabilistic_signal_model, 'covariance')
        args{2} = probabilistic_signal_model.covariance;
        if isfield(probabilistic_signal_model, 'expectation')
            args{3} = probabilistic_signal_model.expectation;
        end
    end
end

[estimate, error_estimate] = do_fuzzy_reduction(measurement, ...
        operator, ideal_operator, args);
    
    function [estimate, error_estimate] = do_fuzzy_reduction(measurement, ...
        operator, ideal_operator, args)
    
    numvarargs = length(args);
    if numvarargs > 3
        error('libfuzzy:probabilistic_reduction:too_many_inputs', ...
            'requires at most 3 optional inputs')
        
    elseif numvarargs == 1
        % Compute a more accurate estimate using the noise covariance matrix.
        noise_covariance = args{1};
        
        tolerance = 1e-8;
        [right_singular singular_values left_singular] = svd(operator' * inv(noise_covariance) * operator);
        inverse_singular_values = diag(pinv(singular_values, tolerance));
        max_inverse_value_index = floor(3.0 / 4.0 * length(inverse_singular_values));
        max_inverse_value = inverse_singular_values(max_inverse_value_index);
        truncated_singular_values = inverse_singular_values <= max_inverse_value;
        singular_values = diag(truncated_singular_values .* inverse_singular_values);
        truncated_operator_inverse = right_singular * singular_values * right_singular';
        reduction_operator = ideal_operator * truncated_operator_inverse * operator';

        estimate = reduction_operator * measurement;
        error_estimate = trace(operator' * inv(noise_covariance) * operator);

    else 
        % Compute basic least squares estimate.
        estimate = pinv(operator) * measurement;
        error_estimate = trace(operator' * operator);
    end

    end

end

