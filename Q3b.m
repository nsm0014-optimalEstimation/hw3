%% Formatting
clc
clear
close all
format shortg

%% Question 3 Part A

% Reading data, defining constants, initializing filter
data = importdata('data/hw3_3.txt');

time = data.data(:,1);
East = data.data(:,2);
North = data.data(:,3);
Psi = data.data(:,4);
Gyro = data.data(:,5);
Radar = data.data(:,6);
dt = 1/5;

H = [diag([ones(1,3)]) zeros(3,2)];

Qd = diag([0.01 0.01 0.01 0 0]);

R = Qd(1:3,1:3)*0.01;

X_mea = [East, North, Psi, zeros(length(time),2)]';
Y_mea = H*X_mea;

X = zeros(5, length(data));
U = data.data(:,5:6)';

X(:,1) = [East(1), North(1), Psi(1), 0 0]';

P_p = zeros(5,5);


for i = 2:length(time)

    Ad = [1   0   0     -sin(X(3,i-1))*dt      0;...
          0   1   0     -cos(X(3,i-1))*dt      0;...
          0   0   1            0             -dt;...
          0   0   0	         1               0;...
          0   0   0            0               1];

    Bd = [sin(X(3,i-1))*dt  0; ...
          cos(X(3,i-1))*dt  0; ...
          0     dt; ...
          0         0; ...
          0         0];


    x_m = Ad*X(:,i-1) + Bd*U(:,i-1);
    P_m = Ad*P_p*Ad' + Qd;

    L = P_m*H'*(H*P_m*H' + R)^-1;
    P_p = (eye(length(x_m)) - L*H)*P_m;
    X(:,i) = x_m + L*(Y_mea(:,i) - H*x_m);

end

Y_est = H*X;

%% Least Squares

H_ls = [1 0 0 0 0 0 0; ...
        0 1 0 0 0 0 0; ...
        0 0 1 0 0 0 0; ...
        0 0 0 1 0 1 0; ...
        0 0 0 0 1 0 1];

R_ls = diag(10*ones(5,1));

P_k = diag(ones(7,1)).^2;

x_k = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
params = zeros(7, length(data));

params(:,1) = x_k;

y = data.data(:,2:6)';

for m = 1:length(time)
    P_k1 = inv(inv(P_k) + H_ls'*inv(R_ls)*H_ls);
    K_k1 = P_k*H_ls'*inv(H_ls*P_k*H_ls' + R_ls);
    x_k1 = x_k + K_k1*(y(:,m) - H_ls*x_k);

    params(:,m) = x_k1;
    P_k = P_k1; K_k = K_k1; x_k = x_k1;
end

fig = figure('Position',[500 250 1000 1000]);
tiledlayout(5,1)
nexttile
hold on
plot(time,params(1,:),'LineWidth',1.5,'LineStyle','--','Color','#398538')
plot(time,X(1,:),'LineWidth',2,'Color','#1b79df')
plot(time,East,'LineWidth',1.5,'LineStyle','-.','Color','#df961b')
ylabel('Position [m]','FontSize',14)
title('East','FontSize',14)
legend('Least Squares','Estimate','Truth',Location='bestoutside',fontsize=14)
xlim([0 60])

nexttile
hold on
plot(time,params(2,:),'LineWidth',1.5,'LineStyle','--','Color','#398538')
plot(time,X(2,:),'LineWidth',2,'Color','#1b79df')
plot(time,North,'LineWidth',1.5,'LineStyle','-.','Color','#df961b')
ylabel('Position [m]','FontSize',14)
title('North','FontSize',14)
xlim([0 60])

nexttile
hold on
plot(time,params(3,:).*180/pi,'LineWidth',1.5,'LineStyle','--','Color','#398538')
plot(time,X(3,:).*180/pi,'LineWidth',2,'Color','#1b79df')
plot(time,Psi*180/pi,'LineWidth',1.5,'LineStyle','-.','Color','#df961b')
ylabel('Rotation [deg]','FontSize',14)
title('Heading','FontSize',14)
xlim([0 60])
    
nexttile
hold on
plot(time,params(4,:),'LineWidth',1.5,'LineStyle','--','Color','#398538')
plot(time,X(4,:).*180/pi,'LineWidth',2,'Color','#1b79df')
ylabel('Bias [m/s]','FontSize',14)
title('Radar bias','FontSize',14)
xlim([0 60])

nexttile
hold on
plot(time,params(5,:),'LineWidth',1.5,'LineStyle','--','Color','#398538')
plot(time,X(5,:),'LineWidth',2,'Color','#1b79df')
ylabel('Bias [m/s]','FontSize',14)
title('Gyroscope bias','FontSize',14)
xlabel('Time [s]','FontSize',14)
xlim([0 60])

saveas(fig,'Q3filter_b.png')