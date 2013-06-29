function [measurement, signal] = probabilistic_measure(operator, ...
                signal_model, noise_model)
%PROBABILISTIC_MEASURE Compute measurement for probabilistic model
%   Detailed explanation goes here

signal = probabilistic_random_vector(signal_model);
noise = probabilistic_random_vector(noise_model);
measurement = operator * signal + noise;

end

