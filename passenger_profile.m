function [values] = passenger_profile(num_values1,mean_value1,std_deviation1,num_values2,mean_value2,std_deviation2)    
    % Generate normally distributed random values
    values1 = mean_value1 + std_deviation1 * randn(1, num_values1);
    values2 = mean_value2 + std_deviation2 * randn(1, num_values2);
    
    values = [values1 values2];
    values = values(values >= 5.5 & values<=23);  
    % Plot the histogram
    figure;
    histogram(values, 'Normalization', 'count');
    title('Passenger profile during the day');
    xlabel('Hour');
    ylabel('Number of Passengers');
end