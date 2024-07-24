% Template for Kalman Filter (KF), Extended Kalman Filter (EKF) and Iterated Extended Kalman Filter (IEKF)

clear; clear all; close all;

% Read data files
load('data_AircraftClimb.mat')
z_k=y_k;

% Initialization
Ts=dt;     %time step (already provided by the data)
N=length(t);      %total number of steps
stdw=[1 5*pi/180];     %standard deviation of w
stdv=[10 0.05];      %standard deviation of v
stdx_0=[20 2 2*pi/180 10];  %standard deviation of x_0
% stdx_0=[200 200 200*pi/180 10];  %standard deviation of x_0
Ex_0=[0 120 0 19000];      %expected value of x_0

xhat_km1_km1 = Ex_0; % x(0|0) = E{x_0}
P_km1_km1 = diag(stdx_0.^2);  % P(0|0) = P(0)
Q=diag(stdw.^2);
R=diag(stdv.^2);

n = length(xhat_km1_km1); % n: state dimension
m = size(u_k, 2);     % m: observation dimension
p = size(z_k, 2);     % m: observation dimension
u_km1 = [zeros(1,m); u_k]; % shifted to have the right indices

% Preallocate storage
stdx_cor  = zeros(N, n);  % \sigma(k-1|k-1), standard deviation of state estimation error (hint: diagonal elements of P(k-1|k-1))
x_cor     = zeros(N, n);  % \hat{x}(k-1|k-1), previous estimation
K       = cell(N, 1);   % K(k) Kalman Gain
innov     = zeros(N, p);  % y(k)-y(k|k-1), innovation, with y(k|k-1)=h(\hat{x}(k|k-1),u(k|k-1),k);

for k=1:N
    % Step 1: Prediction
    [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
    xhat_k_km1 = x_nonlin(end,:); % x(k|k-1) (prediction)

    % Step 2: Covariance matrix of state prediction error / Minimum prediction MSE
    [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts); % Phi(k,k-1), Gamma(k,k-1)
    P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1'; % P(k|k-1) (prediction)

    % Step 3: Kalman Gain
    H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
    Ve = (H_k * P_k_km1 * H_k' + R); % Pz(k|k-1) (prediction)
    K_k = P_k_km1 * H_k' / Ve; % K(k) (gain)
    
    % Step 4: Measurement Update (Correction)
    z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]); % z(k|k-1) (prediction)
    xhat_k_k = xhat_k_km1 + (z_k(k,:) - z_k_km1)*K_k'; % x(k|k) (correction)

    % Step 5: Correction for Covariance matrix of state Estimate error /
    % Minimum MSE
    I_KH = eye(n) - K_k * H_k;
    P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k'; % P(k|k) (correction)

    % Save data: State estimate and std dev
    stdx_cor(k,:) = sqrt(diag(P_km1_km1)); % \sigma(k-1|k-1) Standard deviation of state estimation error
    x_cor(k,:) = xhat_km1_km1; % \hat{x}(k-1|k-1), estimated state
    K{k,1} = K_k; % K(k) (gain)
    innov(k,:)= z_k(k,:) - z_k_km1;
    
    % Recursive step
    xhat_km1_km1 = xhat_k_k; 
    P_km1_km1 = P_k_k;
end


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


%% Functions

function x_dot = funcf(x,u,t)
    
    % Extract values
    h = x(1) ;
    v = x(2) ;
    gamma = x(3);
    m = x(4);
    u1 = u(1) ;
    u2 = u(2) ;
    omega_u1 = 0;
    omega_u2 = 0;
    
    h_dot = v*sin(gamma);
    v_dot = ((- (135*u2^2)/2 - 3/4)*v^2 + u1*cos(u2))/m - (15944000000000*sin(gamma))/1627229345641;
    gamma_dot = cos(gamma)*(v/6378145 - 15944000000000/(1627229345641*v)) + ((175*u2*v^2)/2 + u1*sin(u2))/(m*v);
    m_dot = -u1/15691;

    x_dot = [h_dot; v_dot; gamma_dot; m_dot];

end

function y = funch(x,u,t)

    % Extract values
    h = x(1) ;
    v = x(2) ;
    gamma = x(3);
    m = x(4);
    u1 = u(1) ;
    u2 = u(2) ;
    omega_u1 = 0;
    omega_u2 = 0;
    
    y = [v tan(gamma)];

end

function [Phi,Gamma] = funcLinDisDyn(x_linpt,u_linpt,Ts)

    % Extract values
    h = x_linpt(1) ;
    v = x_linpt(2) ;
    gamma = x_linpt(3);
    m = x_linpt(4);
    u1_m = u_linpt(1) ;
    u2_m = u_linpt(2) ;
    omega_u1 = 0;
    omega_u2 = 0;
        
    % Numerical evaluation of continuous - time dynamics
    A = [[0,                                                                                                               sin(gamma),                                               v*cos(gamma),                                                    0]
        [0,                                                                                          -(2*v*((135*u2_m^2)/2 + 3/4))/m,                 -(15944000000000*cos(gamma))/1627229345641, -((- (135*u2_m^2)/2 - 3/4)*v^2 + u1_m*cos(u2_m))/m^2]
        [0, (175*u2_m)/m + cos(gamma)*(15944000000000/(1627229345641*v^2) + 1/6378145) - ((175*u2_m*v^2)/2 + u1_m*sin(u2_m))/(m*v^2), -sin(gamma)*(v/6378145 - 15944000000000/(1627229345641*v)),         -((175*u2_m*v^2)/2 + u1_m*sin(u2_m))/(m^2*v)]
        [0,                                                                                                                        0,                                                          0,                                                    0]];
    G = [[                         0,                                                                          0]
        [   -cos(omega_u2 - u2_m)/m, ((135*u2_m - 135*omega_u2)*v^2 + sin(omega_u2 - u2_m)*(omega_u1 - u1_m))/m]
        [sin(omega_u2 - u2_m)/(m*v),               (cos(omega_u2 - u2_m)*(omega_u1 - u1_m) - (175*v^2)/2)/(m*v)]
        [                   1/15691,                                                                          0]];

    
    % Discretisation of dynamics
    [Phi , Gamma ]= c2d (A ,G , Ts ) ;

end

function H = funcLinDisObs(x_linpt,u_linpt,t)

    % Extract values
    h = x_linpt(1) ;
    v = x_linpt(2) ;
    gamma = x_linpt(3);
    m = x_linpt(4);
    u1_m = u_linpt(1) ;
    u2_m = u_linpt(2) ;
    omega_u1 = 0;
    omega_u2 = 0;

    % Numerical evaluation
    H = [[0, 1,                0, 0]
        [0, 0, tan(gamma)^2 + 1, 0]];

end
