% MC_FOS_QUALITY - Monte-Carlo simulation to determine quality of
% first order sensor.

MC_ITERATIONS = 1e5;
DIMENSIONS    = 4;

MIN_ALPHA    = -2;
MAX_ALPHA    = +2;
ALPHA_VALUES =  10;

MIN_BETA     = -2;
MAX_BETA     = +2;
BETA_VALUES  =  10;

alphas = linspace(MIN_ALPHA, MAX_ALPHA, ALPHA_VALUES);
betas = linspace(MIN_BETA, MAX_BETA, BETA_VALUES);

%SIGNAL_VARIANCE = sqrt(10.0);
%signal_expectation = sine_waves(DIMENSIONS);
y = 1. / DIMENSIONS : 1. / DIMENSIONS : 1;
signal_expectation = ensure_column(sqrt(2) .* sin(pi .* (1 - y)./2));
%signal_cov_matrix = SIGNAL_VARIANCE * eye(DIMENSIONS);

%probabilistic_signal_model = probabilistic_model(signal_expectation, signal_cov_matrix);
probabilistic_signal_model = probabilistic_model(signal_expectation);
fuzzy_signal_model = coordinated_fuzzy_model(probabilistic_signal_model);

NOISE_STD = 10.0;
noise_expectation = zeros(DIMENSIONS, 1);
noise_cov_matrix = NOISE_STD * eye(DIMENSIONS);
probabilistic_noise_model = probabilistic_model(noise_expectation, noise_cov_matrix);
fuzzy_noise_model = coordinated_fuzzy_model(probabilistic_noise_model);

errors_by_params_probability = zeros(ALPHA_VALUES, BETA_VALUES, MC_ITERATIONS);
errors_by_params_fuzzy = zeros(ALPHA_VALUES, BETA_VALUES, MC_ITERATIONS);

ideal_operator = projection_operator(signal_expectation) ./ DIMENSIONS;

for alpha_index = 1:length(alphas)
    alpha = alphas(alpha_index);
    for beta_index = 1:length(betas)
        beta = betas(beta_index);
        operator = first_order_sensor(alpha, beta, DIMENSIONS);
        for iteration = 1:MC_ITERATIONS
            [probabilistic_measurement, probabilistic_signal_realization] = ...
                probabilistic_measure(operator, probabilistic_signal_model, probabilistic_noise_model);
            [estimate, error_estimate] = probabilistic_reduction(probabilistic_measurement, ...
                operator, ideal_operator, probabilistic_signal_model, probabilistic_noise_model);
            error = (norm(estimate - ideal_operator * probabilistic_signal_realization))^2 / DIMENSIONS;
            errors_by_params_probability(alpha_index, beta_index, iteration) = error;
            
            [fuzzy_measurement, fuzzy_signal_realization] = ...
                fuzzy_measure(operator, fuzzy_signal_model, fuzzy_noise_model);
            [estimate, error_estimate] = fuzzy_reduction(fuzzy_measurement, operator, ...
                ideal_operator, fuzzy_signal_model, fuzzy_noise_model);
            error = (norm(estimate - ideal_operator * fuzzy_signal_realization))^2 / DIMENSIONS;
            errors_by_params_fuzzy(alpha_index, beta_index, iteration) = error;
        end
    end
end

mean_errors_by_params_probability = mean(errors_by_params_probability, 3);
mean_errors_by_params_fuzzy = mean(errors_by_params_fuzzy, 3);

std_errors_by_params_probability = std(errors_by_params_probability, 1, 3);
std_errors_by_params_fuzzy = std(errors_by_params_fuzzy, 1, 3);

figure(1);
surf(betas, alphas, mean_errors_by_params_probability);
%colormap(gray);

figure(2);
surf(betas, alphas, mean_errors_by_params_fuzzy);
%colormap(gray);
