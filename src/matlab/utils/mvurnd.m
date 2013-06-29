function vector = mvurnd(expectation, covariance)
%MVURND Uniformly distributed random numbers
%   MVURND generates a vector of uniformly distributed random numbers with
%   given expected values and covariance matrix.

% iidrv = rand(size(expectation));
% vector = iidrv * (covariance**(1/2));

dimensions = length(expectation);
transformed_covariance = 2 * sin((pi / 6) * covariance);
vector = expectation + normcdf(randn(1, dimensions) * chol(transformed_covariance));
  
end

