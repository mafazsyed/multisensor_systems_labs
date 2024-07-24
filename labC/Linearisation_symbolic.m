clear variables

%% Substitution of the input measurement model.

% define symbolic variables
syms h v gamma m u1 u2 u1_m omega_u1 u2_m omega_u2

% original continuous-time nonlinear system equations
rho = 1;
D = (50/2)*(0.03 + 2.7*u2^2)*rho*v^2;
L = (50/2)*3.5*u2*rho*v^2;

h_dot = v.*sin(gamma);
v_dot = (u1.*cos(u2)-D)./m - 3.9860e+14.*sin(gamma)./6378145.^2;
gamma_dot = (u1.*sin(u2)+L)./(m.*v) + cos(gamma).*(v./6378145-3.9860e+14/(v.*6378145.^2));;
m_dot = -u1./1.5691e+05;

% Final dynamics
x_dot = [h_dot v_dot gamma_dot m_dot];
y = [v tan(gamma)];

% % substitute the input measurement equation in the system equations
x_dot = subs(x_dot,[u1 u2],[u1_m-omega_u1 u2_m-omega_u2]);


% Linearisation of the system 

% find the part corresponding to f
f_expression = subs (x_dot,[omega_u1 omega_u2],[0 0])

% Jacobian linearization to obtain the expression for A
A_expression = jacobian (f_expression, [h v gamma m])

% Jacobian linearization to obtain the expression for B
B_expression = jacobian (f_expression,[u1_m u2_m])

% Jacobian linearization to obtain the expression for G
G_expression = jacobian (x_dot,[omega_u1 omega_u2])

% Jacobian linearization to obtain the expression for H
H_expression = jacobian (y, [h v gamma m])
