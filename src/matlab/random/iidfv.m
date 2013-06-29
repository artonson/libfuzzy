function variables = iidfv(probabilities, count)

%IIDFV Generate independent fuzzy variables
%   IIDFV generates count independent identically distributed fuzzy
%   variables according to probabilities distribution. count should be a
%   positive integer, and probabilities is an array if size 1 ot more. The
%   variables returned are independent in possibilistic sense, and
%   maximally coordinated with probability distribution encoded in
%   probabilities variable.

if count <= 0
    error('count should be given and nonnegative')
end
if length(probabilities) <= 0
    error('probabilities should be nonempty')
end

probabilities = sort(probabilities, 'descend');

number_of_events = length(probabilities);
number_of_vars = count;

variables = zeros(number_of_vars, 1);

% Implementing the naive algorithm for calculation of marginal and 
% conditional distributions for variables, taking O(m^N * N^2) time, 
% where m is the number of simple events, N the number of variables.

for current_variable = 1:number_of_vars
    distribution = zeros(1, number_of_events);
    for current_event = 1:number_of_events
        stack = ones(number_of_vars, 1);
        stack(1:current_variable - 1) = variables(1:current_variable - 1);
        stack(current_variable) = current_event;
        finished_current_event = false;
        while ~finished_current_event
            while stack(number_of_vars) < number_of_events + 1
                distribution(current_event) = ...
                    distribution(current_event) + min(probabilities(stack));
                stack(number_of_vars) = stack(number_of_vars) + 1;
            end
            variable_index = number_of_vars;
            if variable_index == current_variable
                finished_current_event = true;
            else 
                while stack(variable_index) > number_of_events
                    stack(variable_index) = 1;
                    variable_index = variable_index - 1;
                    if variable_index == current_variable
                        finished_current_event = true;
                        break
                    end
                    stack(variable_index) = stack(variable_index) + 1;
                end
            end
        end
    end

    distribution = distribution / sum(distribution);
    variables(current_variable) = find(mnrnd(1, distribution) ~= 0);
end
  
end

