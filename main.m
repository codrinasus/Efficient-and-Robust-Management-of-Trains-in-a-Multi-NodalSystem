M = 1.73e5; %Masa tren[tone]
C = 984; %Capacitate pasageri tren
delta = 30;%Pas esantionare[s]
T = 0:delta:63000;%Orizont timp[s]
I = 1:30;%Index statii metrou
K = 15;%Numar de trenuri in circulatie
di_min = 30;%timp imbarcare minim[s]
di_max = 90;%timp imbarcare maxim[s]
bi_max = 60;%Viteza maxima de imbarcare[pasageri/min]
v_c = 60;%Viteza tren[km/h]
%Coeficienti Davis pentru fortele de frecare
theta_a = -0.0001;
theba_b = -0.25;
theta_c = 30;
%%

% Define time intervals and stations and boarding capacity
delta           = 60;%Sampling time[s]
start_interval  = 7;%Hour at which we want to start the observations
end_interval    = 7.5;%Hour at which we want end the observations
time_horizon    = (end_interval-start_interval)*3600;%Time Horizon
T               = 0:delta:time_horizon; % Set of time intervals 
stations_number = 30; % Number of stations
K               = 18; %Number of trains in service
I               = 1:stations_number; % Set of stations
bi_max         = 120; % Maximum allowable boarding passengers[passengers/min]
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
x  = binvar(stations_number,stations_number,length(T),length(T),K,'full');%space time arc for train presence
bi = intvar(stations_number,length(T),'full');%Boarding passengers at given station
ni = intvar(stations_number,length(T),'full');%Waiting passengers at given station

Constraints = [];
for i = 1 : stations_number
    for t = 2 : length(T)
        Constraints = [Constraints, ni(i,t)>=0];
        Constraints = [Constraints, bi(i,t)>=0];
        Constraints = [Constraints, sum(bi(i,2:t),'all') <= sum(di(i,2:t),'all')];
        Constraints = [Constraints, bi(i,t) <= sum((x(i,i,t-1,t,:).*bi_max),'all')];
    end
end

%Waiting time Objective function

WT = 0;

for i = 1 : stations_number
    for t = 2 : length(T)
        WT = WT + delta * sum(di(1:t,i)-bi(1:t,i),'all');
    end
end

ops = sdpsettings('verbose',2,'debug',0);
optimize(Constraints,WT,ops)
