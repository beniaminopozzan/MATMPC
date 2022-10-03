function xNext = irisStateTransitionFcn(states, control)

Ts = 10e-3;

g = 9.81; % gravity [m/s^2]

mR = 1.5; % drone mass [Kg]
JR = 0.15; % drone inertia
c_t_raw = 5.84e-06 + (-0.1e-6); % thrust coefficient [N/s^2]
c_tau_raw = c_t_raw * 0.06; % drag coefficient [Nm/s^2]
rotor_max_vel = 1100; % [rad/s]
c_direction = [-1 -1 1 1]; % 1: CW , -1: CCW
p_mot = ... % rotor position in body frame [m]
    [0.13 -0.22 0;
    -0.13 0.2 0;
    0.13 0.22 0;
    -0.13 -0.2 0];
theta = atan2(p_mot(:,2),p_mot(:,1)); % arm direction angle
l = vecnorm(p_mot,2,2); % arm length

c_t = c_t_raw .* rotor_max_vel.^2;
c_tau = c_tau_raw .* rotor_max_vel.^2;

p = states(1:3);
q = states(4:7);
v = states(8:10);
omega = states(11:13);
u = states(14:17); % propeller angular rate square
du = control(1:4);

% quaternion conversions
eta = q(1);
epsilon = q(2:4);
epsilon_x = [0, -epsilon(3), epsilon(2); epsilon(3), 0, -epsilon(1); -epsilon(2), epsilon(1), 0]; 
Rq = eye(3) + 2*eta*epsilon_x + 2*epsilon_x^2;
Mq = [eta, -epsilon'; epsilon, eta*eye(3)+epsilon_x];

% thrust matrix

F_temp = [...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    1, 1, 1, 1];
F = 1/mR .* c_t .* F_temp;

% torque matrix
M = c_tau .* F_temp * diag(c_direction) + ...
    c_t .* [...
    sin(theta(1)), sin(theta(2)), sin(theta(3)), sin(theta(4));...
    -cos(theta(1)), -cos(theta(2)), -cos(theta(3)), -cos(theta(4));...
    0, 0, 0, 0]...
    *diag(l);


% explicit ODE RHS
x_dot = [
    v
    1/2*Mq*[0;omega]
    Rq*F*u-[0;0;g]
    1/JR.* M*u
    du
    ];

xNext = states + Ts*x_dot;

end