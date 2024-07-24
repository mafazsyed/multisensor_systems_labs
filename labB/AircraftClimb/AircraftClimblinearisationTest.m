clear variables

% testing values
h=1;
v=140;
gamma=0.01;
m=19000; 
u1=1.5e05;
u2=0.09;

%% Paste your expressions here
A_expression = [[0,                                                                                                       sin(gamma),                                               v*cos(gamma),                                              0]
[0,                                                                                    -(2*v*((135*u2^2)/2 + 3/4))/m,                 -(15944000000000*cos(gamma))/1627229345641, -((- (135*u2^2)/2 - 3/4)*v^2 + u1*cos(u2))/m^2]
[0, (175*u2)/m + cos(gamma)*(15944000000000/(1627229345641*v^2) + 1/6378145) - ((175*u2*v^2)/2 + u1*sin(u2))/(m*v^2), -sin(gamma)*(v/6378145 - 15944000000000/(1627229345641*v)),         -((175*u2*v^2)/2 + u1*sin(u2))/(m^2*v)]
[0,                                                                                                                0,                                                          0,                                              0]
 ];

B_expression = [[            0,                                0]
[    cos(u2)/m,     -(135*u2*v^2 + u1*sin(u2))/m]
[sin(u2)/(m*v), ((175*v^2)/2 + u1*cos(u2))/(m*v)]
[     -1/15691,                                0]];

C_expression = [[0, 1,                0, 0]
[0, 0, tan(gamma)^2 + 1, 0]];

%% Call the test script
AircraftClimbTestReference(h, v, gamma, m, u1, u2 , A_expression, B_expression, C_expression);