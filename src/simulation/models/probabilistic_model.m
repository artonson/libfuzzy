function model = probabilistic_model(expectation, varargin)
%PROBABILISTIC_MODEL Construct probabilistic model of random vector in R^n
%   PROBABILISTIC_MODEL creates a model given variable-length argument list
%   given in the function parameters. With one parameter, it creates a
%   model of a nonrandom vector. With two arguments, it creates a model of
%   a random normally distributed vector. With three argument, a model of
%   random vector from a given distribution is created (currently 'normal' 
%   is supported).

numvarargs = length(varargin);

if nargin < 1 || nargin > 3
    error('libfuzzy:probabilistic_model:too_many_inputs', ...
        'requires at most 3 optional inputs')

elseif numvarargs == 1
    % Create normally distributed vector using expected value and
    % covariance matrix supplied in varargin
    covariance_matrix = varargin{1};
    model = struct('expectation', {expectation}, ...
        'covariance', {covariance_matrix}, ...
        'distribution', {'normal'});

elseif numvarargs == 2
    % Also include the distribution type 
    covariance_matrix = varargin{1};
    distribution = varargin{2};
    model = struct('expectation', {expectation}, ...
        'covariance', {covariance_matrix}, ...
        'distribution', {distribution});

else
    % Create basic model using just expectation
    model = struct('expectation', {expectation});
    
end

end
