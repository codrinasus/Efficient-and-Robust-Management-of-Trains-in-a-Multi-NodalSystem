clear all
close all
clc
%yalmip('clear')
%% Subway Network Parameters
delta             = 60;%Sampling time[s]
start_interval    = 7;%Hour at which we want to start the observations
end_interval      = 9;%Hour at which we want end the observations
K                 = 8; %Number of trains in service
I                 = 17; % Stations number
bi_max            = 120; % Maximum allowable boarding passengers[passengers/min]
%% Train Parameters
Capacity          = 984; %Train capacity
M                 = 1.73e5; %Train mass[t]
v_cruise          = 60; %Subway cruise speed[km/h]
acc               = 1.25; %Subway acceleration[m/s^2]
brk               = 1.2; %Subway breaking[m/s^2]
station_distances = [1.24, 1.65, 1.27, 1.92, 1.55, 1.3, 1.74, 1.13, 1.26, 1.81, 1.87, 1.06, 1.5,1.29]; %List of distance between each station[km]
theta_a           = 0.000136;%Davis parameter for friction
theta_b           = 0.0145;%Davis parameter for friction
theta_c           = 1.244;%Davis parameter for friction
di_min            = 1;%Minimum dwelling time
di_max            = 2;%Maximum dwelling time
TA                = 3;%Time needed at turnaround station 3*delta
omega_e           = 0.5;%Energy cost weight
omega_t           = 1-omega_e;%Waiting time weight
 %% Call train schedule
for i = 1 : length(omega_t)
    [Z(i), WT(i), E(i), x_wait, TD, TE, bi, T, di, Description] = train_schedule(delta, start_interval, end_interval, ...
                                                        K ,I, bi_max, Capacity, M, v_cruise, acc, brk, station_distances,...
                                                        theta_a, theta_b, theta_c, di_min, di_max, TA, omega_t(i), omega_e(i));
end
%% Plotting
%% Passengers boarding
figure()
stairs(sum(value(bi)),'LineWidth',3)
title('Passengers boarding on Subway Line')
xlabel('Time interval [minutes]')
ylabel('Number of passengers')

figure()
hold on
legend_labels = cell(1, K);
%% Train trajectories
figure()
hold on
for k = 1 : K
    trajectory = zeros(length(T),1);
    for t = 1 : length(T)
        for i = 2 : I-1
            trajectory(t) = trajectory(t) + double((i-1)*value(squeeze(x_wait(i,k,t))));
        end
    end

    time = 1 : length(T);
    indices = trajectory >= 1;
    time = time(indices);
    trajectory = trajectory(indices);
    plot(time,trajectory, 'LineWidth', 3)
    legend_labels{k} = ['Train ' num2str(k)];
end
title('Train trajectories')
ylabel('Station number')
xlabel('Time[minutes]')
legend(legend_labels);
hold off
%% Time between train arrivals
i = 5;
Arrivals = zeros(K-1,1);
figure
for k = 2 : K
    Arrivals(k-1) = (value(TE(i,k)) - value(TE(i,k-1)));
end
stairs(Arrivals,'Marker','square','LineWidth',2)
title('Time between train arrivals at station 6')
xlabel('Arrival window')
ylabel('Time between arrivals [minutes]')
%% Savin workspace
save('Simulation.mat')
% %% Function cost versus weights
% figure()
% hold on
% plot(omega_t,value(Z),'LineWidth',3)
% plot(omega_t,value(WT).*omega_t,'LineWidth',3)
% plot(omega_t,value(E).*omega_e,'LineWidth',3)
% title('\omega_T influence in the total and separate costs')
% xlabel('\omega_T')
% ylabel('Cost value')
% legend('Total Cost','Waiting Time','Energy Cost')
% %% X_wait example:
% figure()
% subplot(3,1,1)
% stem(squeeze(value(x_enter(3,1,1:20))),'LineWidth',3)
% xlabel('Time')
% title('x_{enter}')
% 
% subplot(3,1,2)
% stem(squeeze(value(x_exit(3,1,1:20))),'LineWidth',3)
% xlabel('Time')
% title('x_{exit}')
% 
% subplot(3,1,3)
% stem(squeeze(value(x_wait(3,1,1:20))),'LineWidth',3)
% xlabel('Time')
% title('x_{wait}')
% %% Passengers boarding vs demand
% i = 7;
% figure
% stairs(value(bi(i,:)))
% hold on
% stairs(value(di(i,:)))
% title('Time between train arrivals at station 6')
% xlabel('Arrival window')
% ylabel('Time between arrivals [minutes]')