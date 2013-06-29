function fos = first_order_sensor(alpha, beta, dimensions)
%FIRST_ORDER_SENSOR Matrix for first order sensor
%   FIRST_ORDER_SENSOR creates matrix for the first order sensor given
%   its parameters and dimensions

integration_step = 1 / dimensions;
grid = linspace(0, 1, dimensions + 1);
x1 = grid(2 : dimensions + 1);
x2 = grid(1 : dimensions);

temp = ones(1, dimensions);
t1 = (x1' * temp - temp' * x1) .* beta ./ alpha;
t2 = (x2' * temp - temp' * x1) .* beta ./ alpha;
t3 = (x1' * temp - temp' * x2) .* beta ./ alpha;
t4 = (x2' * temp - temp' * x2) .* beta ./ alpha;

CONST = -alpha / beta^2 / integration_step;
d = CONST.*(exp(-t1)-exp(-t2)-exp(-t3)+exp(-t4));
di = 1 / beta + CONST * (1 - exp(-beta / alpha * integration_step));
fos = d - triu(d) + di .* eye(size(d));

end
