function [estimate, error_estimate] = fuzzy_reduction(measurement, ...
        operator, signal_model, noise_model)
%FUZZY_REDUCTION Summary of this function goes here
%   Detailed explanation goes here

B = Sigma^(-1/2) * operator * F^(1/2);
y = Sigma^(-1/2) * (measurement - operator * myf0);

exit_flat = 0;
iterations = 0;
MAX_ITERATIONS = 10;
max_value = 2;
while exit_flat ~= 1 && iterations <= MAX_ITERATIONS
    [w, value, exit_flat] = fzero(@(w) fuzzy_reduction_helper(w,B,y), [0 max_value]);
    iterations = iterations + 1;
end





R = F * operator' * inv(operator * F * operator' + w * Sigma);
R_xi = f0 + R * (xi - operator * f0);
h = trace(F - R * operator * F);

    function f = fuzzy_reduction_helper(x, B, y)
    b = B'*y;
    A = B' * B + x * eye(size(B' * B));
    z2 = A \ b;

    n1 = norm(y - B*z2);
    n2 = norm(z2);
    f = n1 - n2;

    % T = B * B';
    % f = x * norm(     inv(T + x * eye(size(T))) \ * y) - ...
    %         norm(B' * inv(T + x * eye(size(T))) * y);
    % 

    end

end

