function vector = fuzzy_random_vector(model)
%FUZZY_RANDOM_VECTOR Fuzzy random vector from model
%   FUZZY_RANDOM_VECTOR creates a random vector according to the
%   supplied fuzzy model. (See src/simulation/models/fuzzy_model.m 
%   for reference concerning avaiable model parameters.)
%   Number of granules is assumed to be the same over all events.

if isfield(model, 'coordinates')
    number_of_vars = length(model.coordinates);
    number_of_events = length(model.coordinates(1).granules);
    probabilities = prepare_probabilities(ones(number_of_events, 1), number_of_vars);
    vars = iidfv(probabilities, number_of_vars);
    vector = zeros(number_of_vars, 1);
    for coordinate = 1:number_of_vars
        vector(coordinate) = random_point_from_granules(...
            model.coordinates(coordinate).granules, vars(coordinate));
    end
else
    vector = model.expectation;
end

end

    function point = random_point_from_granules(granules, index)
    intervals = granules(index).intervals;
    [rows, cols] = size(intervals);
    equal_probabilities = ones(rows, 1) / rows;
    interval_index = find(mnrnd(1, equal_probabilities) ~= 0);

    pos_inf_index = find(intervals == +Inf);
    neg_inf_index = find(intervals == -Inf);
    if ~isempty(pos_inf_index)
        intervals(pos_inf_index) = +999999;
    end
    if ~isempty(neg_inf_index)
        intervals(neg_inf_index) = -999999;
    end
    point = intervals(interval_index, 1) + ...
        (intervals(interval_index, 2) - intervals(interval_index, 1)) * rand;
    end
