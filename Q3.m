%% Formatting
clc
clear
close all
format shortg

%% Question 3 Part A

% Reading data, defining constants, initializing filter
data = importdata("data\hw3_3.txt");

time = data.data(:,1);
eastPos = data.data(:,2);
northPos = data.data(:,3);
radarBias = data.data(:,4);
psi = data.data(:,5);
gyroBias = data.data(:,6);

A_c = 

Y = [eastPos northPos radarBias psi gyroBias];
C = [1 1 0 1 0];

dt = diff(time);
dt = dt(1);

num = length(time);





