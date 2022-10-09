%% Formatting
clc
clear
close all
format shortg
%% Question 2 Part A

% Reading data, defining constants, initializing filter
data = importdata("data\hw3_2.txt");

plotFontSize = 14;

time = data(:,1);
dt = diff(time);
Y_data = data(:,2);

A_c = 0;
Q = 0;
B = 0;
C = 1;
trick = expm([-A_c (B*Q*B');zeros(1) A_c']*dt(1));
A_d = trick(2,2);
Q_d = A_d*trick(1,2);
R_d = 1;

% Pre-allocating for speed
L_k = zeros(1,length(time));
X_hat = zeros(1,length(time));
X = 0;
P_plus = 1;
% Beginning Kalman Filter Simulation
for i = 1:length(time)
    % Time Update
    X = A_d*X;
    P_minus = A_d*P_plus*A_d' + Q_d;

    % Measurement Update
    L_k(i) = P_minus*C'*(C*P_minus*C' + R_d)^-1;
    P_plus  = (1 - L_k(i)*C)*P_minus;
    X_hat(i) = X + L_k(i)*(Y_data(i) - C*X);


    X = X_hat(i);
end

P_plus_ss = ((A_d*P_plus*A_d' + Q_d)^-1 + C*R_d^-1*C')^-1
L_k_ss = P_plus_ss*C'*R_d^-1

% Plotting Simulation
simFig = figure('Position',[500 250 1000 600],'Name','Estimator vs Data simulation');
hold on
plot(time,Y_data,'LineWidth',2)
plot(time,X_hat,'LineWidth',2)
legend('Measurement','Filtered Position')
ylabel('Meters [m]')
xlabel('Time [s]')
fontsize(simFig,plotFontSize,"points")
saveas(simFig,'Q2_simulation.png')

% Plotting Simulation
KFFig = figure('Position',[500 250 1000 600],'Name','Estimator vs Data simulation');
hold on
plot(time,L_k,'LineWidth',2)
ylabel('Gain Value')
xlabel('Time [s]')
fontsize(KFFig,plotFontSize,"points")
saveas(KFFig,'Q2_gains.png')

%% Question 2 Part B
clc
clear

% Reading data, defining constants, initializing filter
data = importdata("data\hw3_2.txt");

plotFontSize = 14;

time = data(:,1);
dt = diff(time);
Y_data = data(:,2);

A_c = 0;
Q = [0.0001 0.001 0.01 0.1];

% Pre-allocating for speed
L_k = zeros(4,length(time));
X_hat = zeros(4,length(time));

% Beginning Kalman Filter Simulation
for j = 1:length(Q)
    X = 0;
    P_plus = 1;
    B = 1;
    C = 1;
    trick = expm([-A_c (B*Q(j)*B');zeros(1) A_c']*dt(1));
    A_d = trick(2,2);
    Q_d = A_d*trick(1,2);
    R_d = 1;


    for i = 1:length(time)
        % Time Update
        X = A_d*X;
        P_minus = A_d*P_plus*A_d' + Q_d;

        % Measurement Update
        L_k(j,i) = P_minus*C'*(C*P_minus*C' + R_d)^-1;
        P_plus  = (1 - L_k(j,i)*C)*P_minus;
        X_hat(j,i) = X + L_k(j,i)*(Y_data(i) - C*X);


        X = X_hat(j,i);
    end
    P_plus_ss(j) = ((A_d*P_plus*A_d' + Q_d)^-1 + C*R_d^-1*C')^-1;
    L_k_ss(j) = P_plus_ss(j)*C'*R_d^-1;
end


% Plotting Simulation
simBFig = figure('Position',[500 250 1000 600],'Name','Estimator vs Data simulation');
hold on
plot(time,Y_data,'LineWidth',2)
plot(time,X_hat,'LineWidth',2)
legend('Measurement','Q = 0.0001','Q = 0.001','Q = 0.01','Q = 0.1','Location','best')
ylabel('Meters [m]')
xlabel('Time [s]')
fontsize(simBFig,plotFontSize,"points")
saveas(simBFig,'Q2b_simulation.png')

% Plotting Simulation
KFBFig = figure('Position',[500 250 1000 600],'Name','Estimator vs Data simulation');
hold on
plot(time,L_k,'LineWidth',2)
ylabel('Gain Value')
xlabel('Time [s]')
legend('Q = 0.0001','Q = 0.001','Q = 0.01','Q = 0.1','Location','best')
fontsize(KFBFig,plotFontSize,"points")
saveas(KFBFig,'Q2b_gains.png')

%% Problem 2 Part C

yf = filter(sqrt(Q_d),[1 -(1-sqrt(Q_d))],Y_data,Y_data(1));

% Plotting Simulation
simCFig = figure('Position',[500 250 1000 600],'Name','Estimator vs Data simulation');
hold on
plot(time,Y_data,'LineWidth',2)
plot(time,X_hat(4,:),'LineWidth',2)
plot(time,yf,'LineWidth',1,'LineStyle','--')
legend('Measurement','Kalman Filter Position','>>filter position','Location','best')
ylabel('Meters [m]')
xlabel('Time [s]')
fontsize(simCFig,plotFontSize,"points")
saveas(simCFig,'Q2c_simulation.png')