
s = sprintf('events\tvars\ttime');
disp(s);
for number_of_events = 1:10
    for number_of_vars = 1:10
        tic;
        probabilities = prepare_probabilities(ones(number_of_events,1), number_of_vars);
        vars = iidfv(probabilities, number_of_vars);
        time_to_generate = toc;
        s = sprintf('%d\t%d\t%g', number_of_events, number_of_vars, time_to_generate);
        disp(s);
    end
end