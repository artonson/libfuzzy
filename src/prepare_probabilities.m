function probabilities = prepare_probabilities(ordering, number_of_variables)
%PREPARE_PROBABILITIES Compute probability distribution for iidfv function
%   Function IIDFV requires a valid probability distribution to function.
%   As of now (15 Jun 2012 - @artonson), IIDFV can generate fuzzy variables
%   from a possibility distribution in which ordering of the possibilities
%   is 0.1..1, i.e. p[i] > p[i+1] for any i. 
%   To stochastically model such a possibility distribution (for example, 
%   to generate N independent fuzzy variables from it), one needs a
%   special probability distribution. In it, probability of an event 
%   (w[i[1]], ..., w[i[N]]) is the value of the minimum of 
%   probabilities of its constituents w[i[k]].
%
%   PREPARE_PROBABILITIES gives out a probability distribution described in
%   [Pytyev Zhivotnikov morphology paper].

number_of_events = length(ordering);
probabilities = zeros(number_of_events, 1);

numerator_value = 1;
for current_event = 1:number_of_events - 1
    if current_event > 1
        numerator_value = numerator_value - ((current_event - 1) ^ number_of_variables - ...
            (current_event - 2) ^ number_of_variables) * probabilities(current_event - 1);
    end
    denominator_value = current_event ^ number_of_variables - (current_event - 1) ^ number_of_variables;
    lhs_value = numerator_value / (denominator_value + 1);
    rhs_value = numerator_value / denominator_value;
    probabilities(current_event) = lhs_value + (rhs_value - lhs_value) / 2.;
end
if number_of_events > 1
	numerator_value = numerator_value - ((number_of_events - 1) ^ number_of_variables - ...
        (number_of_events - 2) ^ number_of_variables) * probabilities(number_of_events - 1);
end
probabilities(number_of_events) = numerator_value / (number_of_events ^ number_of_variables - (number_of_events - 1) ^ number_of_variables);

end
