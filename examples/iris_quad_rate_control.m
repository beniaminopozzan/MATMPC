%------------------------------------------%

% coplanar quadrotor based on IRIS model

%------------------------------------------%


%% Dimensions

nx=13;  % No. of differential states
nu=4;  % No. of controls
nz=0;  % No. of algebraic states
ny=17; % No. of outputs
nyN=13; % No. of outputs at the terminal point
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

rollTimeConstant = 0.0729; % [s]
pitchTimeConstant = 0.09068; % [s]
yawTimeConstant = 0.1515; % [s]

p=states(1:3);
q=states(4:7);
v=states(8:10);
omega=states(11:13);

u=controls(1:4); % thrust & roll-pitch-yaw rates
thrust = u(1);
roll_rate_setpoint = u(2);
pitch_rate_setpoint = u(3);
yaw_rate_setpoint = u(4);

% quaternion conversions
eta = q(1);
epsilon = q(2:4);
epsilon_x = [
    0, -epsilon(3), epsilon(2);
    epsilon(3), 0, -epsilon(1);
    -epsilon(2), epsilon(1), 0]; 
Rq = eye(3) + 2*eta*epsilon_x + 2*epsilon_x^2;
Mq = [eta, -epsilon';
    epsilon, eta*eye(3)+epsilon_x];

% explicit ODE RHS
x_dot = [
    v
    1/2*Mq*[0;omega]
    Rq./mR*[0; 0; thrust] - [0;0;g]
    1./rollTimeConstant*(roll_rate_setpoint-omega(1))
    1./pitchTimeConstant*(pitch_rate_setpoint-omega(2))
    1./yawTimeConstant*(yaw_rate_setpoint-omega(3))
    ];

% algebraic function
z_fun = [];

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;

% integratori diversi (vedi orbitali)
     
%% Objectives and constraints

roll = 0;%atan2(2*(q(1)*q(2)+q(3)*q(4)),1-2*(q(2)^2+q(3)^2));
pitch = 0;%asin(2*(q(1)*q(3)-q(4)*q(2)));
yaw = 0;%atan2(2*(q(1)*q(4)+q(2)*q(3)),1-2*(q(3)^2+q(4)^2));

% inner objectives
h = [
    p
    roll
    pitch
    yaw
    q'*refsq
    v
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
