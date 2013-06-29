function signal = sine_waves(dimensions)
%SINE_WAVES Generate two-component sinusoidal signal
%   SINE_WAVES generates a signal consisting of two components where the
%   first component is a large-magnitude low-frequency signal, and the
%   second if low-magintude high-frequency signal. The first one
%   corresponds to the overall shape of the signal, whereas the second one
%   describes `thin` signal details.

large_magnidute = 1;
small_magnidute = 0.1 * large_magnidute;

low_frequency = 3 * pi / dimensions;
high_frequency = 10 * low_frequency;

signal = (large_magnidute * sin((1:dimensions) * low_frequency) + ...
            small_magnidute * sin((1:dimensions) * high_frequency))';

end

