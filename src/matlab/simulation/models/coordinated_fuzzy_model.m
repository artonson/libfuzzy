function fuzzy_model = coordinated_fuzzy_model(probabilistic_model)
%COORDINATED_FUZZY_MODEL Create a fuzzy model coordinated with given
%probabilistic model
%   COORDINATED_FUZZY_MODEL creates a fuzzy model maximally coordinated
%   with a given probabilistic model.

sharpness = 1.1;
ordering = 1:3;

if sharpness <= 1.0
    error('libfuzzy:coordinated_fuzzy_model:small_sharpness', ...
        'sharpness parameter should be greater then 1.0')
end

fuzzy_model.coordinated_probabilistic_model = probabilistic_model;

if isfield(probabilistic_model, 'distribution')
    switch probabilistic_model.distribution
        case 'normal'
            variances = diag(probabilistic_model.covariance);
            for coordinate = 1:length(probabilistic_model.expectation)
                expectation = probabilistic_model.expectation(coordinate);
                variance = variances(coordinate);
                fuzzy_model.coordinates(coordinate).granules = ...
                    granulate_normal_pdf(expectation, variance, ordering, sharpness);
            end
        otherwise
            error('libfuzzy:coordinated_fuzzy_model:unknown_distribution', ...
                'unknown distribution')
    end
else
    fuzzy_model.expectation = probabilistic_model.expectation;
end

end
