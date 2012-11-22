function measurement = fuzzy_measure(operator, signal_model, noise_model)
%MEASURE_FUZZY Compute measurement for fuzzy model
%   Detailed explanation goes here

signal = fuzzy_random_vector(signal_model);
noise = fuzzy_random_vector(noise_model);
measurement = operator * signal + noise;

end

