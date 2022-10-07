%% Formatting
clc
clear
close all
format shortg
%% Begin Problem 1 Part A
% Constants/Simulation Parameters
tFinal = 100; % [s]
tStart = 0; % [s]
dt = 0.1; % [s]

tSequence = tStart:dt:tFinal; % [s] (array)

plotFontSize = 14;

% State Space System
A_CL = [0 1;-1 -1.4]; % Closed-Loop A-Matrix
A_D = expm(A_CL); % X = e^(At)
B = [0;1];
C = [1 0];

% Pre-allocating for speed
X = zeros(2,length(tSequence));
Y = zeros(1,length(tSequence));
w = zeros(1,length(tSequence));
Q = zeros(1,length(tSequence));
eta = zeros(1,length(tSequence));
R = zeros(1,length(tSequence));

% Beginning Simulation
for i = 1:length(tSequence)
    
    % Simulating Process Noise
    w(i) = 2*randn();
    Q(i) = mean(w(i)*w(i)');

    % Simulating Sensor Noise
    eta(i) = randn();
    R(i) = mean(eta(i)*eta(i)'); 
    % Simulation is Continuous
    X_dot = A_CL*X(:,i) + B*w(i);
    X(:,i+1) = X_dot*dt + X(:,i);

    % Measurement
    Y(i+1) = C*X(:,i+1) + randn(1,1);
end

% Plotting Simulation
simFig = figure('Position',[500 500 1000 600],'Name','Continuous Model simulation');
hold on
plot(tSequence,Y(2:end),'LineWidth',2)
plot(tSequence,X(1,2:end),'LineWidth',2)
legend('Measured Position','Actual Position')
ylabel('Meters [m]')
xlabel('Time [s]')
fontsize(simFig,plotFontSize,"points")
saveas(simFig,'Q1_simulation.png')

%% Begin Problem 1 Part B
Q = mean(Q);
Q_d = B*Q*B'*dt;
R_d = mean(R);
