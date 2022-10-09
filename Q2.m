%% Formatting 
clc
clear
close all
format shortg
%% Question 2 Part A

% Reading data, defining constants, initializing filter
data = importdata("data\hw3_2.txt");

time = data(:,1);
dt = diff(time);
Y = data(:,2);

