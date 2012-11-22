% MC_FOS_QUALITY - Monte-Carlo simulation to determine quality of
% first order sensor.

MC_ITERATIONS =  1000;
DIMENSIONS    = 4;

MIN_ALPHA    = -3;
MAX_ALPHA    = +3;
ALPHA_VALUES =  20;

MIN_BETA     = -3;
MAX_BETA     = +3;
BETA_VALUES  =  20;

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

NOISE_VARIANCE = 10.0;
noise_expectation = zeros(DIMENSIONS, 1);
noise_cov_matrix = NOISE_VARIANCE * eye(DIMENSIONS);
probabilistic_noise_model = probabilistic_model(noise_expectation, noise_cov_matrix);
fuzzy_noise_model = coordinated_fuzzy_model(probabilistic_noise_model);

errors_by_params_probability = zeros(ALPHA_VALUES, BETA_VALUES);
errors_by_params_fuzzy = zeros(ALPHA_VALUES, BETA_VALUES);

ideal_operator = projection_operator(signal_expectation) ./ DIMENSIONS;

for alpha_index = 1:length(alphas)
    alpha = alphas(alpha_index);
    for beta_index = 1:length(betas)
        beta = betas(beta_index);
        operator = first_order_sensor(alpha, beta, DIMENSIONS);
        for iteration = 1:MC_ITERATIONS
%             [probabilistic_measurement, probabilistic_signal_realization] = ...
%                 probabilistic_measure(operator, probabilistic_signal_model, probabilistic_noise_model);
%             [estimate, error_estimate] = probabilistic_reduction(probabilistic_measurement, ...
%                 operator, ideal_operator, probabilistic_signal_model, probabilistic_noise_model);
%             error = (norm(estimate - ideal_operator * probabilistic_signal_realization))^2 / DIMENSIONS;
%             errors_by_params_probability(alpha_index, beta_index) = ...
%                 errors_by_params_probability(alpha_index, beta_index) + error;
            tic;
            [fuzzy_measurement, fuzzy_signal_realization] = ...
                fuzzy_measure(operator, fuzzy_signal_model, fuzzy_noise_model);
            [estimate, error_estimate] = fuzzy_reduction(fuzzy_measurement, operator, ...
                ideal_operator, fuzzy_signal_model, fuzzy_noise_model);
            error = (norm(estimate - ideal_operator * fuzzy_signal_realization))^2 / DIMENSIONS;
            if isnan(error)
                disp('est'); disp(estimate);
                disp(fuzzy_signal_realization);
                disp(fuzzy_measurement);
            end
            errors_by_params_fuzzy(alpha_index, beta_index) = ...
                errors_by_params_fuzzy(alpha_index, beta_index) + error;
            toc;
        end
    end
end
% errors_by_params_probability = errors_by_params_probability / MC_ITERATIONS;
errors_by_params_fuzzy = errors_by_params_fuzzy / MC_ITERATIONS;

% figure(1);
% surfl(betas, alphas, errors_by_params_probability);
% shading interp;
% colormap(pink);

figure(2);
surfl(betas, alphas, errors_by_params_fuzzy);
shading interp;
colormap(pink);

