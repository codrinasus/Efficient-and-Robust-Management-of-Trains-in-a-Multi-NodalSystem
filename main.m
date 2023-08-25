%% Waiting Time

% Define time intervals and stations and boarding capacity
delta           = 30;%Sampling time[s]
start_interval  = 7;%Hour at which we want to start the observations
end_interval    = 7.25 ;%Hour at which we want end the observations
time_horizon    = (end_interval-start_interval)*3600;%Time Horizon
T               = 0:delta:time_horizon; % Set of time intervals 
stations_number = 30; % Number of stations
K               = 6; %Number of trains in service
I               = 1:stations_number; % Set of stations
bi_max          = 120; % Maximum allowable boarding passengers[passengers/min]
tolerance       = 1e-16;%Constraints tolerance

%Dynamic passenger demand
p_ij_c  = assign_passenger_demand(start_interval, end_interval, delta);
p_ij    = zeros(stations_number,stations_number,length(T));
for t=1:length(T)
    for i = 1 : stations_number
        for j = 1 : stations_number
            if i < j && ...
               ((i <= stations_number/2) && (j<= stations_number/2) || ...
               (i>= stations_number/2) && (j>= stations_number/2))
                p_ij(j,i,t) = p_ij_c(t);
            else
                p_ij(j,i,t) = 0;
            end
        end
    end
end
  
%Time variant passengers demand
di = zeros(stations_number, length(T));

for t=1:length(T)
    for i = 1 : stations_number
        for j = 1 : stations_number
            di(i,t) = di(i,t) + p_ij(j,i,t) * delta;
        end
    end
end

% Create YALMIP variables
x  = binvar(stations_number,stations_number,length(T),K,'full');%space time arc for train presence
bi = intvar(stations_number,length(T),'full');%Boarding passengers at given station

%Waiting time constraints (28) && (29) && (31)
Constraints = [];
for i = 1 : stations_number
    for t = 1 : length(T)
        Constraints = [Constraints, bi(i,t)>=0];
        Constraints = [Constraints, sum(bi(i,1:t)) <= sum(di(i,1:t))];
        Constraints = [Constraints, bi(i,t) <= sum((x(i,i,t,:).*bi_max))];
    end
end

%Waiting time Objective function

WT = 0;

for i = 1 : stations_number
    for t = 1 : length(T)
        WT = WT + delta * sum(di(i,1:t)-bi(i,1:t));
    end
end

%% Energy consumption
Capacity          = 984; %Train capacity
M                 = 1.73e5; %Train mass[t]
v_cruise          = 60; %Subway cruise speed[km/h]
acc               = 1.25; %Subway acceleration[m/s^2]
brk               = 1.2; %Subway breaking[m/s^2]
station_distances = [1.13,1.5,1.16,1.84,1.41,1.19,1.59,1.03,1.15,1.65,1.7,0.97,1.37,1.18]; %List of distance between each station[km]
acc_time          = v_cruise * 1000 / (3600 * acc); %Time needed for train to reach cruise speed[s]
brk_time          = v_cruise * 1000 / (3600 * brk); %Time needed for train to break[s]
theta_a           = -0.0001;%Davis parameter for friction
theta_b           = -0.25;%Davis parameter for friction
theta_c           = 30;%Davis parameter for friction
di_min            = 30;%Minimum dwelling time
di_max            = 90;%Maximum dwelling time
T_a               = acc_time;%acceleration timestamp
T_b               = brk_time;%breaking timestamp
T_c               = zeros(stations_number-2,1);%cruising timestamp
for i = 1 : stations_number/2-1
    T_c(i) = (station_distances(i)/v_cruise)*3600;
end
T_c(stations_number/2:end) = T_c(stations_number/2-1:-1:1);

%Dwelling time bounds constraints (24) && (25)
for i = 2 : stations_number-1
    for k = 1 : K
        for t = 1 : length(T)-di_max/delta-1
              Constraints = [Constraints, x(i-1,i,t,k) >= x(i,i+1,(t+1):(t+di_min/delta),k)];
              Constraints = [Constraints, x(i-1,i,t,k) <= x(i,i+1,t+1+di_max/delta,k)];
        end
    end
end

h = ceil((T_c+acc_time+brk_time)/delta);%headway time on a link
%Time distance between adjacent trains  (22) && (23) 
for i = 1 : stations_number - 1
    if(i == 15)
                continue
    end
    for k = 2 : K
        for t = 1 : length(T) - max(h)
            if(i < 15)
                Constraints = [Constraints, x(i,i+1,(t+1):(t+h(i)),k-1)<= x(i,i+1,t,k)];
            else
                Constraints = [Constraints, x(i,i+1,(t+1):(t+h(i-1)),k-1)<= x(i,i+1,t,k)];
            end
        end
    end
end


TR = 6;%Time to prepare for next cycle 6*delta
TA = 3;%Time needed at turnaround station 3*delta

%Cycle constraint (21)
for k = 1 : K
    for t = 1 : length(T) - TA
        Constraints = [Constraints, x(15,15,(t+1:t+TA),k) <= x(15,15,t,k)];
    end
end 

%Turnaround constraint (20)
for k = 1 : K
    for t = 1 : length(T) - TR
        Constraints = [Constraints, x(15,15,(t+1:t+TR),k) <= x(30,30,t,k)];
    end
end

%Concurency station
for i = 1 : stations_number
        for t = 1 : length(T)
            Constraints = [Constraints, sum(x(i,i,t,:)) == 1];
        end
end

%Concurency link
for i = 1 : stations_number - 1
    for j = i + 1 : stations_number
        for t = 1 : length(T)
            Constraints = [Constraints, sum(x(i,j,t,:)) == 1];
        end
    end
end

% %Present only in one place at given time
% for t = 1 : length(T)
%     for k = 1 : K
%         Constraints = [Constraints, sum(x(:,:,t,k),'all') == 1];
%     end
% end

%Passenger capacity constraint (30)
for i = 1 : stations_number 
    for t = 1 : length(T)
        insum = 0;
        for k = 1 : K
            insum = insum + sum(di(i,t)*(1-x(i,i,t,k)));
        end
            Constraints = [Constraints, insum <= Capacity];
    end
end

E = 0;%Energy cost (17)

for i = 1 : stations_number -1
    if i == 15
        continue
    end
    for j = i : stations_number
        for t = 1 : length(T)
            for k = 1 : k
                E = E + M * (acc + 1/M * (theta_a * v_cruise)^2 + theta_b *v_cruise + theta_c) * ...
                    v_cruise * delta * x(i,j,t,k);
            end
        end
    end
end

omega_t = 0.5;%Waiting time weight
omega_e = 1 - omega_t;%Energy cost weight

Z = omega_t * WT + omega_e * E;%objective function (18)

%Solving the problem
ops = sdpsettings('verbose',2,'debug',1);
optimize(Constraints,Z,ops);

% Plotting
figure()
stairs(sum(value(bi)))
title('Passengers boarding on Subway Line')
xlabel('Time interval [minutes]')
ylabel('Number of passengers')

figure()
hold on
for i = 1 : stations_number - 1
    plot(2*i + squeeze(value(x(i,i+1,:,6))))
end
hold off