%% Formatting
clc
clear
close all
format shortg

%% Begin Problem 4

% Defining parameters given in problem statement
A_c = [-2.62 12;-0.96 -2];
B_c = [14;1];
C_c = [1 0];



Q_c = 0.05;

dt = 0.1; % 10 [Hz]
t = 0:dt:5;

measurement_noise = gaussianDistFCN([length(t) 1],0.1,0);
% Using Bryson's Trick to convert the continuous system to discrete
BT = [-A_c B_c*Q_c*B_c'; zeros(2) A_c'];
BT = expm(BT*dt);

A_da = BT(3:4,3:4)';
Q_d = A_da*BT(1:2,3:4);
% B_d = B_c*dt;
C_d = C_c;
R_d = 0.05;

steer_angle = [zeros(5,1); ones(46,1)];

% Starting the simulation
x = zeros(2,length(t));
y = zeros(1,length(t));

x_plus = zeros(2,length(t));
P_plus = zeros(2,2);

% Part b - updated dynamic model
A_c = [-2.42 4;-0.99 -2];
B_c = [18;1];
C_c = [1 0];




BT = [-A_c B_c*Q_c*B_c'; zeros(2) A_c'];
BT = expm(BT*dt);

A_d = BT(3:4,3:4)';
B_d = B_c*dt;

for i = 2:length(t)

    % Modeling System
    x(:,i) = A_d*x(:,i-1) + B_d*steer_angle(i-1);
    y(i) = C_d*x(:,i) + measurement_noise(i);

    % Modeling Observer
    x_minus = A_da*x_plus(:,i-1);
    P_minus = A_da*P_plus*A_da' + Q_d;

    L = P_minus*C_d'*(C_d*P_minus*C_d' + R_d)^-1;
    P_plus = (eye(2) - L*C_d)*P_minus;
    x_plus(:,i) = x_minus + L*(y(i) - C_d*x_minus);


end
L_ss_poles = eig(A_c - L*C_d)
%% Plotting
figQ4b = figure('Position',[250 100 1000 600]);
tiledlayout(2,1)

nexttile
hold on
plot(t,x_plus(1,:),'LineWidth',2)
plot(t,x(1,:),'LineWidth',2)
ylabel('Yaw Rate')
legend('Estimate','Truth',Location='southeast')
ax = gca;
ax.FontSize = 16;


nexttile
hold on
plot(t,x_plus(2,:),'LineWidth',2)
plot(t,x(2,:),'LineWidth',2)
ylabel('Sideslip')
xlabel('Time [s]')
legend('Estimate','Truth',Location='northeast')
ax = gca;
ax.FontSize = 16;

saveas(figQ4b,'Q4b.png')
