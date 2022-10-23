%% Formatting
clc
clear
close all
format shortg

%% Question 3 Part A

% Reading data, defining constants, initializing filter
data = importdata("data\hw3_3.txt");

%% Optimal Estimation - Homework 3 - Problem 3

clear 
clc
close all

%% Part A

data = importdata('data/hw3_3.txt');

% Time Initialization
time = data.data(:,1);
k = length(data.data(:,1));

% Saving data to variables
East = data.data(:,2);
North = data.data(:,3);
psi = data.data(:,4);
psi_dot = data.data(:,5);
velo = data.data(:,6);

% EKF Initialization
xhat = zeros(7,k);
xhat(1,1) = East(1);
xhat(2,1) = North(1);
xhat(3,1) = psi(1);
xhat(4,1) = psi_dot(1);
xhat(5,1) = velo(1);

procNoise = 0.01;
measNoise = .01;

P = eye(7); 

Q = procNoise * eye(7);
Q(6,6) = 0.00001;
Q(7,7) = 0.00001;

R = measNoise^2 * eye(5);

H = [eye(5) zeros(5,2)];
H(4,6) = 1;
H(5,7) = 1;

P_log = cell(k,1);
P_log{1} = P;

% EKF
for i = 1:k-1

    % Time Step Calculation
    dt = time(i+1) - time(i);

    % Time Update
    xhat(1,i+1) = xhat(1,i) + xhat(5,i)*sin(xhat(3,i))*dt;
    xhat(2,i+1) = xhat(2,i) + xhat(5,i)*cos(xhat(3,i))*dt;
    xhat(3,i+1) = xhat(3,i) + xhat(4,i)*dt;
    xhat(4,i+1) = xhat(4,i) - xhat(6,i);
    xhat(5,i+1) = xhat(5,i) - xhat(7,i);
    xhat(6,i+1) = xhat(6,i);
    xhat(7,i+1) = xhat(7,i);

    A = eye(7);
    A(1,3) = xhat(5,i)*cos(xhat(3,i))*dt;
    A(1,5) = sin(xhat(3,i))*dt;
    A(2,3) = -xhat(5,i)*sin(xhat(3,i))*dt;
    A(2,5) = cos(xhat(3,i))*dt;
    A(3,4) = dt;
    A(4,6) = -1;
    A(5,7) = -1;

    P = A*P*A' + Q;
    
    % Measurement Update
    
    if i >= 1300
        H = zeros(5,7);
    end
    
    K = P*H'*(H*P*H' + R)^-1;
    
    z = [East(i+1); North(i+1); psi(i+1); psi_dot(i+1);
        velo(i+1)];

    xhat(:,i+1) = xhat(:,i+1) + K*(z - H*xhat(:,i+1));

    P = (eye(7) - K*H)*P;

    P_log{i+1} = P;

end

fig = figure('Position',[250 0 1000 1000]);
tiledlayout(5,1)
nexttile
hold on
plot(time,xhat(2,:),'LineWidth',2)
plot(time,East,'LineWidth',2)
ylabel('Position [m]','FontSize',14)
title('Estimated East','FontSize',14)

nexttile
hold on
plot(time,xhat(3,:),'LineWidth',2)
plot(time,North,'LineWidth',2)
ylabel('Position [m]','FontSize',14)
title('Estimated North','FontSize',14)

nexttile
hold on
plot(time,xhat(5,:).*180/pi,'LineWidth',2)
plot(time,psi,'LineWidth',2)
ylabel('Rotation [deg]','FontSize',14)
title('Estimated Yaw','FontSize',14)

nexttile
hold on
plot(time,xhat(4,:).*180/pi,'LineWidth',2)
plot(time,psi_dot,'LineWidth',2)
ylabel('Rotation Rate [deg/s]','FontSize',14)
title('Estimated Yaw Rate','FontSize',14)

nexttile
hold on
plot(time,xhat(6,:),'LineWidth',2)
plot(time,velo,'LineWidth',2)
ylabel('Velocity [m/s]','FontSize',14)
xlabel('Time [s]','FontSize',14)
title('Estimated Velocity','FontSize',14)



figure
plot(time,East,'.')
hold on
plot(time,xhat(1,:),'.')
title('East (m)')

figure
plot(time,North,'.')
hold on
plot(time,xhat(2,:),'.')
title('North (m)')

figure
plot(time,psi,'.')
hold on
plot(time,xhat(3,:),'.')
title('Heading (rad)')

figure
plot(time,psi_dot,'.')
hold on
plot(time,xhat(4,:),'.')
title('Gyro (rad/s)')

figure
plot(time,velo,'.')
hold on
plot(time,xhat(5,:),'.')
title('Radar (m/s)')

figure
plot(time,xhat(6,:))
hold on
plot(time,xhat(7,:))
title('Mesaurement Biases')
legend('Gyro Bias','Radar Bias')

%% Part B
a = 0;
b = 0;
for i = 1:199
    dt = time(i+1) - time(i);
    a = a + psi_dot(i+1)*dt;
    b = b + velo(i+1)*dt;
end






