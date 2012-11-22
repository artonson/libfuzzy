function vector = probabilistic_random_vector(model)
%PROBABILISTIC_RANDOM_VECTOR Random vector from model
%   PROBABILISTIC_RANDOM_VECTOR creates a random vector according to the
%   supplied stochastic model. (See 
%   src/simulation/models/probabilistic_model.m for reference concerning 
%   avaiable model parameters.)

if isfield(model, 'covariance')
    if isfield(model, 'distribution')
        if strcmp(model.distribution, 'normal')
            vector = (mvnrnd(model.expectation, model.covariance))';
        elseif strcmp(model.distribution, 'uniform')
            vector = (mvnrnd(model.expectation, model.covariance))';
        else
            error('libfuzzy:probabilistic_random_vector:wrong_distribution', ...
                'given distribution is unknown')
        end
    else
        vector = (mvnrnd(model.expectation, model.covariance))';
    end    
else
    vector = model.expectation;
end
vector = ensure_column(vector);
end

