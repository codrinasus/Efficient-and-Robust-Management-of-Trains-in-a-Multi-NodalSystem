function p_ij_c = assign_passenger_demand(start_interval,end_interval,delta)
    %start_interval and end_interval are both specified in integer
    %format[hours]
    %The time step(delta) is specified in seconds [s]
    passenger_p = passenger_profile(28000,7.5,1.5,65000,16.5,3);%passenger profile
    close()
    passenger_p = passenger_p(passenger_p >= start_interval & passenger_p <= end_interval+delta/3600);%restrain the passenger profile to the specified time interval
    
    k = 1;
    p_ij_c = zeros(1,(end_interval-start_interval)/(delta/3600)+1);
    for t = start_interval:delta/3600:end_interval
        p_ij_c(k) = numel(passenger_p(passenger_p >= t & passenger_p <= t+delta/3600))/(delta*15*15);
        k = k + 1;
    end
end