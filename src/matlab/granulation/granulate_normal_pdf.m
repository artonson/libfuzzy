function granules = granulate_normal_pdf(mean, std, ...
	ordering, sharpness)
%GRANULATE_NORMAL_PDF Create granules for normal density
%   Detailed explanation goes here

if sharpness <= 1.0
    error('libfuzzy:granulate_normal_pdf:small_sharpness', ...
        'sharpness parameter should be greater then 1.0')
end

total_probability = 0;
for granule_index = 1:length(ordering)
    probability = sharpness ^ granule_index / (1 + sharpness) ^ granule_index;
    intervals = zeros(1, 2);
    cdf_probabilities = [(1.0 - probability) / 2.0, (1.0 + probability) / 2.0];
    outer_interval = norminv(cdf_probabilities, mean, std);
    if granule_index == 1
        intervals = outer_interval;
    elseif granule_index == length(ordering)
        intervals(1, :) = [-Inf, granules(granule_index - 1).intervals(1, 1)];
        [rows, cols] = size(granules(granule_index - 1).intervals);
        intervals(2, :) = [granules(granule_index - 1).intervals(rows, 2), +Inf];
        probability = 1 / (1 + sharpness) ^ granule_index;        
    else
        intervals(1, :) = [granules(granule_index - 1).intervals(1, 1), outer_interval(1, 1)];
        [rows, cols] = size(granules(granule_index - 1).intervals);
        intervals(2, :) = [outer_interval(1, 2), granules(granule_index - 1).intervals(rows, 2)];
    end
    granules(granule_index).intervals = intervals;
    granules(granule_index).probability = probability;
    total_probability = total_probability + probability;
end

% normalise
for granule_index = 1:length(ordering)
    granules(granule_index).probability = granules(granule_index).probability / total_probability;
end

end
