function signal = stairs(amplitudes, dimensions)
%STAIRS Generate stairs-like signal with varying amplitudes
%   SINE_WAVES generates a signal consisting of two components where the
%   first component is a large-magnitude low-frequency signal, and the
%   second if low-magintude high-frequency signal. The first one
%   corresponds to the overall shape of the signal, whereas the second one
%   describes `thin` signal details.

if length(amplitudes) >= length(dimensions)
    error('libfuzzy:stairs:too_low_dimensions', ...
        'number of amplitudes should not exceed signal dimensions')
end

signal = ensure_column(zeros(dimensions, 1));




end

