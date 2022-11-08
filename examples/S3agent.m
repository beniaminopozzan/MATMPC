%------------------------------------------%
% Inverted Pendulum 
  
% "Autogenerating microsecond solvers for nonlinear MPC: A tutorial
% using ACADO integrators", Quirynen, 2015

% typical configuration: 1) N=80, Ts=Ts_st=0.025, no shifting 2) N=40,
% Ts=Ts_st=0.05, shifting

%------------------------------------------%


%% Dimensions

nx=7;  % No. of differential states
nu=3;  % No. of controls
nz=0;  % No. of algebraic states
ny=7; % No. of outputs
nyN=4; % No. of outputs at the terminal point
np=4; % No. of model parameters
nc=0; % No. of general constraints
ncN=0; % No. of general constraints at the terminal point
nbx = 0; % No. of bounds on states
nbu = 3; % No. of bounds on controls

% state and control bounds
nbx_idx = []; % indexs of states which are bounded
nbu_idx = 1:3; % indexs of controls which are bounded

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

J = [
    1 0 0
    0 1 0
    0 0 1.5
    ];

q       = states(1:4);
omega   = states(5:7);
u       = controls(1:3);
qref    = params(1:4);

omega_sk = [
    0   -omega(3) omega(2)
    omega(3)    0   -omega(1)
    -omega(2)   omega(1)    0
];
eta = q(1);
epsilon = q(2:4);
epsilon_sk = [
    0   -epsilon(3) epsilon(2)
    epsilon(3)    0   -epsilon(1)
    -epsilon(2)   epsilon(1)    0
];
M = [
    -epsilon.'
    eta*eye(3)+epsilon_sk
    ];
q_dot = 0.5*M*omega;
omega_dot = J\(omega_sk*J*omega+u);

% explicit ODE RHS
x_dot=[q_dot; omega_dot];
 
% algebraic function
z_fun = [];                   

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;        
     
%% Objectives and constraints

% inner objectives

h = [
    acos(q'*qref)
    omega
    u
    ];
hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality constraints
general_con = [];
general_con_N = [];

%% NMPC discretizing time length [s]

Ts_st = 0.01; % shooting interval time
