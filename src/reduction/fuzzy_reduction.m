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
        reduction_operator = ideal_operator * ...
            pinv(operator' * inv(noise_covariance) * operator, 1e-5) *...
            operator';
        estimate = reduction_operator * measurement;
        error_estimate = trace(operator' * inv(noise_covariance) * operator);

    else 
        % Compute basic least squares estimate.
        estimate = pinv(operator) * measurement;
        error_estimate = trace(operator' * operator);
    end

    end

end

