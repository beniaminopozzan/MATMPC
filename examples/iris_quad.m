%------------------------------------------%

% coplanar quadrotor based on IRIS model

%------------------------------------------%


%% Dimensions

nx=13;  % No. of differential states
nu=4;  % No. of controls
nz=0;  % No. of algebraic states
ny=20; % No. of outputs
nyN=16; % No. of outputs at the terminal point
np=0; % No. of model parameters
nc=0;%0; % No. of general inequality constraints
ncN=0;%1; % No. of general inequality constraints
nbx = 0; % No. of bounds on states
nbu = 4; % No. of bounds on controls

% state and control bounds
nbx_idx = [];  % indexs of states which are bounded
nbu_idx = 1:4;  % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);   % differential states
controls = SX.sym('controls',nu,1); % control input
alg      = SX.sym('alg',nz,1);      % algebraic states
params   = SX.sym('paras',np,1);    % parameters
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable


%% Dynamics

%Constant

g = 9.81; % gravity [m/s^2]

mR = 1.5; % drone mass [Kg]
JR = 0.15; % drone inertia
c_t_raw = 5.84e-06; % thrust coefficient [N/s^2]
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

p=states(1:3);
q=states(4:7);
v=states(8:10);
omega=states(11:13);
u=controls(1:4); % propeller angular rate square

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
x_dot = [Rq*v; 1/2*Mq*[0;omega]; F*u-Rq'*[0;0;g]; 1/JR.* M*u];

% algebraic function
z_fun = [];

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

roll = atan2(2*(q(1)*q(2)+q(3)*q(4)),1-2*(q(2)^2+q(3)^2));
pitch = asin(2*(q(1)*q(3)-q(4)*q(2)));
yaw = atan2(2*(q(1)*q(4)+q(2)*q(3)),1-2*(q(3)^2+q(4)^2));

% inner objectives
h = [
    p
    roll
    pitch
    yaw
    q
    v
    omega
    u
    ];
hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality path constraints
general_con = [];
general_con_N = [];

%% NMPC discretizing time length [s]

Ts_st = 0.01; % shooting interval time
