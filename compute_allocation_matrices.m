g = 9.81;

mass = 1.5;

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

% thrust matrix

F_temp = [...
    0, 0, 0, 0;...
    0, 0, 0, 0;...
    1, 1, 1, 1];
F = c_t .* F_temp;

% torque matrix
M = c_tau .* F_temp * diag(c_direction) + ...
    c_t .* [...
    sin(theta(1)), sin(theta(2)), sin(theta(3)), sin(theta(4));...
    -cos(theta(1)), -cos(theta(2)), -cos(theta(3)), -cos(theta(4));...
    0, 0, 0, 0]...
    *diag(l);

allocationMat = [
    F(3,:)
    M
    ];

takeOffThrustCoef = 0.7;
requiredTakeOffThrust = mass*g;

maxThrust = [1 0 0 0]*allocationMat*[1;1;1;1];

allocationMat = allocationMat./maxThrust;

save("data/allocationMat.mat",'allocationMat');