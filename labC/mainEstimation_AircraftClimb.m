% Template for Kalman Filter (KF), Extended Kalman Filter (EKF) and Iterated Extended Kalman Filter (IEKF)

clear; clear all; close all;

% Read data files
load('data_AircraftClimb.mat')
N=length(t);

[x_cor] = runEKF(u_k,y_k,t,dt);

%% Plots

figure
subplot(2,2,1)
plot(t(1:N),actual_state.x_k(1:N,1),'-g',t(1:N),x_cor(1:N,1),'--b','LineWidth',2)
grid on

subplot(2,2,2)
plot(t(1:N),y_k(1:N,1),'-.k',t(1:N),actual_state.x_k(1:N,2),'-g',t(1:N),x_cor(1:N,2),'--b','LineWidth',2)
grid on
legend('raw measurements','original','estimation')

subplot(2,2,3)
plot(t(1:N),actual_state.x_k(1:N,3),'-g',t(1:N),x_cor(1:N,3),'--b','LineWidth',2)
grid on

subplot(2,2,4)
plot(t(1:N),actual_state.x_k(1:N,4),'-g',t(1:N),x_cor(1:N,4),'--b','LineWidth',2)
grid on

