%% Numerical Solution of a ODE Equation - ode45 (Runge-Kutta)
clear all
close all
clc

% Initial conditions
x0 = [0 129 0 19000];

% Inputs (if constant, otherwise specify it in the dynamics function)
u = [1.5e05 0.09];

% Timespan from t0 to tf
t0 = 0;             % initial time
tf = 50;             % final time
tspan = [t0 tf];    % timespan

%% use ode45 for the solution of the ODEs (nonlinear)
[t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u, t), tspan, x0);   %solve the initial value problem
u_nonlin=repmat(u,length(t_nonlin),1);                            %formulate input vector
y_nonlin=funch(x_nonlin, u_nonlin, t_nonlin);                     %formulate output vector

%% Configuration for linearised models and simulations
% Linearisation point
x_linpt=x0;
u_linpt=[0 0];

% Initial output
y0=funch(x_linpt, u_linpt, 0);

%% use ode45 for the solution of the ODEs (linear)

% Solve the initial value problem 
% Note that for the linearised system, the states x and inputs u in x_dot=A*x+B*u represent
% deviation from the linearisation point (x_linpt, u_linpt)
[t_lin, x_lin] = ode45(@(t,x) funcLinDyn(x, transpose(u-u_linpt), t, x_linpt, u_linpt), tspan, x0-x_linpt);
u_lin=repmat(u-u_linpt,length(t_lin),1);                     %Input vector
y_lin=funcLinObs(x_lin, u_lin, t_lin, x_linpt, u_linpt);     %output vector

% Recover actual states, inputs and outputs
x_lin=x_lin+repmat(x_linpt,length(t_lin),1);
u_lin=u_lin+repmat(u_linpt,length(t_lin),1);
y_lin=y_lin+repmat(y0,length(t_lin),1);

%% use loop for the solution of the ODEs (linear discrete time)

% time step
Ts=0.01;
% time vector
t_lindis=t0:Ts:tf;

%state vector (initialization)
x_lindis=zeros(length(t_lindis),length(x0));       
x_lindis(1,:)=x0-x_linpt;

% loop based on x_{k+1}=Phi*x_{k}+Gamma_u*u_{k}
for i=1:length(t_lindis)-1
    x_lindis(i+1,:)=transpose(funcLinDiscDyn(transpose(x_lindis(i,:)), transpose(u-u_linpt), i, x_linpt, u_linpt, Ts));       %state vector (initialization)
end
u_lindis=repmat(u-u_linpt,length(t_lindis),1);       %input vector
y_lindis=funcLinDiscObs(x_lindis, u_lindis, t_lindis, x_linpt, u_linpt);

% Recover actual states, inputs and outputs
x_lindis=x_lindis+repmat(x_linpt,length(t_lindis),1);
u_lindis=u_lindis+repmat(u_linpt,length(t_lindis),1);
y_lindis=y_lindis+repmat(y0,length(t_lindis),1);

%% Plotting the results
figure;
hold on
plot(t_nonlin, x_nonlin(:,1),'LineWidth',2); 
plot(t_lin, x_lin(:,1),'--','LineWidth',2); 
plot(t_lindis, x_lindis(:,1),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('State x1 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

figure;
hold on
plot(t_nonlin, x_nonlin(:,2),'LineWidth',2); 
plot(t_lin, x_lin(:,2),'--','LineWidth',2); 
plot(t_lindis, x_lindis(:,2),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('State x2 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

figure;
hold on
plot(t_nonlin, x_nonlin(:,3),'LineWidth',2); 
plot(t_lin, x_lin(:,3),'--','LineWidth',2); 
plot(t_lindis, x_lindis(:,3),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('State x3 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

figure;
hold on
plot(t_nonlin, x_nonlin(:,4),'LineWidth',2); 
plot(t_lin, x_lin(:,4),'--','LineWidth',2); 
plot(t_lindis, x_lindis(:,4),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('State x4 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on


figure;
hold on
plot(t_nonlin, u_nonlin(:,1),'LineWidth',2); 
plot(t_lin, u_lin(:,1),'--','LineWidth',2); 
plot(t_lindis, u_lindis(:,1),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('Input u1 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

figure;
hold on
plot(t_nonlin, u_nonlin(:,2),'LineWidth',2); 
plot(t_lin, u_lin(:,2),'--','LineWidth',2); 
plot(t_lindis, u_lindis(:,2),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('Input u2 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

figure;
hold on
plot(t_nonlin, y_nonlin(:,1),'LineWidth',2); 
plot(t_lin, y_lin(:,1),'--','LineWidth',2); 
plot(t_lindis, y_lindis(:,1),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('Output y1 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

figure;
hold on
plot(t_nonlin, y_nonlin(:,2),'LineWidth',2); 
plot(t_lin, y_lin(:,2),'--','LineWidth',2); 
plot(t_lindis, y_lindis(:,2),'-.','LineWidth',2); 
xlabel('Time (s)');
ylabel('Output y2 (-)');
legend('Nonlinear Model','Linearised Model','Linearised Discrete-time Model')
grid on

%% Local Functions
% wrapper for the dynamics equation
function x_dot = funcf(x, u, t)
    h=x(1);
    v=x(2);
    gamma=x(3);
    m=x(4);
    u1=u(1);
    u2=u(2);
   
    % Compute derivatives
    rho = 1;
    D = (50/2)*(0.03 + 2.7*u2^2)*rho*v^2;
    L = (50/2)*3.5*u2*rho*v^2;
    
    h_dot = v.*sin(gamma);
    v_dot = (u1.*cos(u2)-D)./m - 3.9860e+14.*sin(gamma)./6378145.^2;
    gamma_dot = (u1.*sin(u2)+L)./(m.*v) + cos(gamma).*(v./6378145-3.9860e+14/(v.*6378145.^2));
    m_dot = -u1./1.5691e+04;
    
    x_dot = [h_dot; v_dot; gamma_dot; m_dot];
end

%wrapper for the measurement equation
function y = funch(x, u, t)
    h=x(:,1);
    v=x(:,2);
    gamma=x(:,3);
    m=x(:,4);
    u1=u(:,1);
    u2=u(:,2);
   
    % Compute output
    y = [v tan(gamma)];
end

%wrapper for the linearised dynamics equation
function x_dot = funcLinDyn(x, u, t, x_linpt, u_linpt)
    %Linearisation point
    h=x_linpt(1);
    v=x_linpt(2);
    gamma=x_linpt(3);
    m=x_linpt(4);
    u1=u_linpt(1);
    u2=u_linpt(2);
   
    % Compute derivatives
    A=[[0,                                                                                                       sin(gamma),                                               v*cos(gamma),                                              0]
    [0,                                                                                    -(2*v*((135*u2^2)/2 + 3/4))/m,                 -(15944000000000*cos(gamma))/1627229345641, -((- (135*u2^2)/2 - 3/4)*v^2 + u1*cos(u2))/m^2]
    [0, (175*u2)/m + cos(gamma)*(15944000000000/(1627229345641*v^2) + 1/6378145) - ((175*u2*v^2)/2 + u1*sin(u2))/(m*v^2), -sin(gamma)*(v/6378145 - 15944000000000/(1627229345641*v)),         -((175*u2*v^2)/2 + u1*sin(u2))/(m^2*v)]
    [0,                                                                                                                0,                                                          0,                                              0]
    ];
    
    B=[[            0,                                0]
    [    cos(u2)/m,     -(135*u2*v^2 + u1*sin(u2))/m]
    [sin(u2)/(m*v), ((175*v^2)/2 + u1*cos(u2))/(m*v)]
    [     -1/15691,                                0]];

    
    x_dot = A*x+B*u;
end

%wrapper for the linearised measurement/observation equation 
function y = funcLinObs(x, u, t, x_linpt, u_linpt)
    h=x_linpt(1);
    v=x_linpt(2);
    gamma=x_linpt(3);
    m=x_linpt(4);
    u1=u_linpt(1);
    u2=u_linpt(2);

    C=[[0, 1,                0, 0]
    [0, 0, tan(gamma)^2 + 1, 0]];
   
    % Compute output
    y = transpose(C*x');
end

%wrapper for the linearised discrete-time dynamics equation
function [x_kp1] = funcLinDiscDyn(x_k, u_k, k, x_linpt, u_linpt, Ts)
    %Linearisation point
    h=x_linpt(1);
    v=x_linpt(2);
    gamma=x_linpt(3);
    m=x_linpt(4);
    u1=u_linpt(1);
    u2=u_linpt(2);
   
    % Compute derivatives
    A=[[0,                                                                                                       sin(gamma),                                               v*cos(gamma),                                              0]
    [0,                                                                                    -(2*v*((135*u2^2)/2 + 3/4))/m,                 -(15944000000000*cos(gamma))/1627229345641, -((- (135*u2^2)/2 - 3/4)*v^2 + u1*cos(u2))/m^2]
    [0, (175*u2)/m + cos(gamma)*(15944000000000/(1627229345641*v^2) + 1/6378145) - ((175*u2*v^2)/2 + u1*sin(u2))/(m*v^2), -sin(gamma)*(v/6378145 - 15944000000000/(1627229345641*v)),         -((175*u2*v^2)/2 + u1*sin(u2))/(m^2*v)]
    [0,                                                                                                                0,                                                          0,                                              0]
    ];
    
    B=[[            0,                                0]
    [    cos(u2)/m,     -(135*u2*v^2 + u1*sin(u2))/m]
    [sin(u2)/(m*v), ((175*v^2)/2 + u1*cos(u2))/(m*v)]
    [     -1/15691,                                0]];


    [Phi,Gamma_u]=c2d(A,B,Ts);
    x_kp1 = Phi*x_k+Gamma_u*u_k;
end

%wrapper for the linearised discrete-time measurement/observation equation 
function y = funcLinDiscObs(x, u, t,x_linpt, u_linpt)
    h=x_linpt(1);
    v=x_linpt(2);
    gamma=x_linpt(3);
    m=x_linpt(4);
    u1=u_linpt(1);
    u2=u_linpt(2);

    H=[[0, 1,                0, 0]
    [0, 0, tan(gamma)^2 + 1, 0]];
   
    % Compute output
    y = transpose(H*x');
end


