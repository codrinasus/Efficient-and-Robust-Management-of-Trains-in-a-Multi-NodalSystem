clear all
close all
clc
yalmip('clear')
%% Subway Network Parameters
delta             = 60;%Sampling time[s]
start_interval    = 7;%Hour at which we want to start the observations
end_interval      = 8;%Hour at which we want end the observations
time_horizon      = (end_interval-start_interval)*3600;%Time Horizon
T                 = 0:delta:time_horizon;%Set of time intervals 
K                 = 6; %Number of trains in service
I                 = 16; % Stations number
bi_max            = 120; % Maximum allowable boarding passengers[passengers/min]
tolerance         = 1e-16;%Constraints tolerance
%% Train Parameters
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
di_min            = 1;%Minimum dwelling time
di_max            = 2;%Maximum dwelling time
T_a               = acc_time;%acceleration timestamp
T_b               = brk_time;%breaking timestamp
T_c               = zeros(I-4,1);%cruising timestamp
TR                = 6;%Time to prepare for next cycle 6*delta
TA                = 3;%Time needed at turnaround station 3*delta
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

%% Travelling time between stations
for i = 1 : (I-4)/2
    T_c(i) = (station_distances(i)/v_cruise)*3600;
end
T_c((I-4)/2+1:end) = T_c((I-4)/2:-1:1);
% h = ceil((T_c+acc_time+brk_time)/delta);%headway time on a link
% h = [h(1:(I-4)/2); TA; h((I-4)/2+1:end)];


%% Create YALMIP variables
x_enter  = binvar(I,K,length(T),'full');%1 if train k enters station i at time t
x_exit   = binvar(I,K,length(T),'full');%1 if train k exits station i at time t
x_wait   = binvar(I,K,length(T),'full');
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
%     Constraints = [Constraints, x_exit(I,k,end) == 1];
    for t = length(T) : -1 : 2
        Constraints = [Constraints, x_exit(I,k,1:t-1) <= x_exit(I,k,t)];
    end
end

%Station is occupied
for i = 2 : I-1
    for k = 1 : K
        Constraints = [Constraints, sum(x_enter(i,k,:)) == 1, ...
                                    sum(x_exit(i,k,:)) == 1 ];
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
        Constraints = [Constraints, TE(i,k) - TE(i,k-1) >= h(i-2)];
    end
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
        Constraints = [Constraints, bi(i-1,t) <= sum((x_wait(i,:,t).*bi_max/delta))];
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
            E = E + M * (acc + 1/M * (theta_a * v_cruise)^2 + theta_b *v_cruise + theta_c) * ...
                v_cruise * delta * x_wait(i,k,t);
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
omega_t = 0.5;%Waiting time weight
omega_e = 1 - omega_t;%Energy cost weight

Z = omega_t * WT + omega_e * E;%objective function (18)

ops = sdpsettings('verbose',2,'debug',1);
optimize(Constraints,Z,ops);

%% Plotting
figure()
stairs(sum(value(bi)))
title('Passengers boarding on Subway Line')
xlabel('Time interval [minutes]')
ylabel('Number of passengers')

for k = 1 : K
    figure()
    hold on
    for i = 2 : I - 1
        plot(i * squeeze(value(x_enter(i,k,:))))
    end
    hold off
end