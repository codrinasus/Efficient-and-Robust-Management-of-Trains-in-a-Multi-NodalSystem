function [Z, WT, E, x_wait, TD, TE, bi,T, di, Description] = train_schedule(delta, start_interval, end_interval, ...
                                                        K ,I, bi_max, Capacity, M, v_cruise, acc, brk, station_distances,...
                                                        theta_a, theta_b, theta_c, di_min, di_max, TA, omega_t, omega_e)
    %% Parameters
    time_horizon      = (end_interval-start_interval)*3600;%Time Horizon
    T                 = 0:delta:time_horizon;%Set of time intervals
    acc_time          = v_cruise * 1000 / (3600 * acc); %Time needed for train to reach cruise speed[s]
    brk_time          = v_cruise * 1000 / (3600 * brk); %Time needed for train to break[s]
    T_a               = acc_time;%acceleration timestamp
    T_b               = brk_time;%breaking timestamp
    T_c               = zeros(I-4,1);%cruising timestamp
    h                 = ones(I-3,1);%Headway time between stations
    %% Dynamic passenger demand
    p_ij_c  = assign_passenger_demand(start_interval, end_interval, delta);
    p_ij    = zeros(I-2,I-2,length(T));
    for t=1:length(T)
        for i = 1 : I-2
            for j = 1 : I-2
                if i < j && ...
                        ((i <= (I-2)/2) && (j<= (I-2)/2) || ...
                        (i>= (I-2)/2) && (j>= (I-2)/2))
                    p_ij(j,i,t) = p_ij_c(t);
                else
                    p_ij(j,i,t) = 0;
                end
            end
        end
    end
    
    %Time variant passengers demand
    di = zeros(I-2, length(T));
    for t=1:length(T)
        for i = 1 : I-2
            for j = 1 : I-2
                di(i,t) = di(i,t) + p_ij(j,i,t) * delta;
            end
        end
    end
    
    di = di.*30;
    %% Travelling time between stations
    for i = 1 : (I-3)
        T_c(i) = (station_distances(i)/v_cruise)*3600;
    end
    %T_c((I-4)/2+1:end) = T_c((I-4)/2:-1:1);
    h = ceil((T_c+acc_time+brk_time)/delta);%headway time on a link
    % h = [h(1:(I-4)/2); TA; h((I-4)/2+1:end)];
    
    
    %% Create decision variables
    x_enter  = binvar(I,K,length(T),'full');%1 if train k enters station i at time t
    x_exit   = binvar(I,K,length(T),'full');%1 if train k exits station i at time t
    x_wait   = binvar(I,K,length(T),'full');%1 if train k waits in station i at time t
    bi       = intvar(I-2,length(T),'full');%Boarding passengers at given station
    TE       = intvar(I,K,'full');%Station entering times
    TD       = intvar(I,K,'full');%Station departure times
    %% Constraints
    
    Constraints = [];
    
    %Dummy start station
    for k = 1 : K
        Constraints = [Constraints, x_enter(1,k,1) == 1];
        for t = 2 : length(T)
            Constraints = [Constraints, x_enter(1,k,(t):end) <= x_enter(1,k,t-1)];
        end
    end
    
    %Dummy stop station
    for k = 1 : K
        %Constraints = [Constraints, x_exit(I,k,end) == 1];
        for t = length(T) : -1 : 2
            Constraints = [Constraints, x_exit(I,k,1:t-1) <= x_exit(I,k,t)];
        end
    end
    
    %Station is occupied
    for i = 2 : I-1
        for k = 1 : K
            Constraints = [Constraints, sum(x_enter(i,k,:)) <= 1, ...
                sum(x_exit(i,k,:)) <= 1 ];
        end
    end
    
    %Train cannot enter a station later than it exits
    for i = 2 : I-1
        for k = 1 : K
            for t = 1 : length(T)
                Constraints = [Constraints, sum(x_enter(i,k,1:t)) >= sum(x_exit(i,k,1:t))];
            end
        end
    end
    
    %Entering Time
    for i = 2 : I-1
        for k = 1 : K
            TE(i,k) = sum(x_enter(i,k,:).* [1:length(T)]);
        end
    end
    
    %Departure Time
    for i = 2 : I-1
        for k = 1 : K
            TD(i,k) = sum(x_exit(i,k,:).* [1:length(T)]);
        end
    end
    
    %Time inbetween succesive stations
    for i = 3 : I-1
        for k = 1 : K
            Constraints = [Constraints, TE(i,k) - TD(i-1,k) >= h(i-2)];
        end
    end
    
    %Time between arrivals and departures
    for i = 2 : I-1
        for k = 1 : K
            Constraints = [Constraints, TD(i,k) - TE(i,k) >= di_min];
            Constraints = [Constraints, TD(i,k) - TE(i,k) <= di_max];
        end
    end
    
    %Time between succesive trains
    for i = 3 : I-1
        for k = 2 : K
            Constraints = [Constraints, TE(i,k) - TD(i,k-1) >= h(i-2)];
        end
    end
    
    %Constraints to let only one train to leave the dummy station
    for k = 2 : K
        Constraints = [Constraints, TE(2,k) - TE(2,1:k-1) >= di_max];
    end
    
    %X_WAIT
    for i = 2 : I - 1
        for k = 1 : K
            for t = 1 : length(T)
                x_wait(i,k,t) = sum(x_enter(i,k,1:t)) - sum(x_exit(i,k,1:t));
            end
        end
    end
    
    
    %Boarding passengers constraints
    for i = 2 : I-1
        for t = 1 : length(T)
            Constraints = [Constraints, bi(i-1,t)>=0];
            Constraints = [Constraints, sum(bi(i-1,1:t)) <= sum(di(i-1,1:t))];
            Constraints = [Constraints, bi(i-1,t) <= sum((x_wait(i,:,t).*bi_max*delta/60))];
        end
    end
    
    %Passenger capacity constraint (30)
    for i = 2 : I-1
        for t = 1 : length(T)
            Constraints = [Constraints, di(i-1,t)*sum((1-x_wait(i,:,t))) <= Capacity];
        end
    end
    
    %% Energy Cost objective function
    E = 0;
    for i = 2 : I -1
        for t = 1 : length(T)
            for k = 1 : k
                E = E + M * (acc + 1/M * (theta_a * v_cruise)^2 + theta_b *v_cruise +...
                    theta_c) * v_cruise * delta * (1-x_wait(i,k,t));
            end
        end
    end
    
    %% Waiting time Objective function
    
    WT = 0;
    for i = 2 : I - 1
        for t = 1 : length(T)
            WT = WT + delta * sum(di(i-1,1:t)-bi(i-1,1:t));
        end
    end
    %% Solving the problem
    Z = omega_t * WT + omega_e * E;%objective function (18)
    
    ops = sdpsettings('verbose',2,'debug',1);
    Description = optimize(Constraints,Z,ops);
end