% General Parameters
addpath('Functions_general');
clc
clear
close all

%%% CHANGE! - DIMENSIONS
nx = 2; %%% NUMBER OF STATES
nu = 1; %%% NUMBER OF INPUTS
nc = 4; %%% NUMBER OF CONSTRAINTS

T = 0.5;
N = 8;

%%% CHANGE! - MATRICES
A = eye(nx);
for i = 1:nx-1
    A(i, i+1) = T;
end

B = zeros(nx, nu);
B(end) = T;

Q = eye(nx);
% A = [1 T 0; 0 1 T; 0 0 1]; %%% CHANGE!
% B = [0; 0; T];             %%% CHANGE!
% Q = [1 0 0; 0 1 0; 0 0 1]; %%% CHANGE!

R = 0.1;        
[K, Px]  = dlqr(A, B, Q, R); 
K        = -K;                     
Phi      = A + B*K; 
% Below: Construct the true uncertainty set of EV and LV
theta_EV = pi/40;           
theta_LV = pi/40;           
Rotation_EV = [cos(theta_EV) -sin(theta_EV); sin(theta_EV) cos(theta_EV)]; 
Rotation_LV = [cos(theta_LV) -sin(theta_LV); sin(theta_LV) cos(theta_LV)]; 
V_EV_rot = (Rotation_EV * [-0.06 -0.015; 0.06 -0.015; 0.01 0.025; -0.01 0.025]')';
V_LV_rot = (Rotation_LV * [-0.06 0.015; 0.06 0.015; 0.01 -0.025; -0.01 -0.025]')';

Xi_true_EV = embedInFullDim(Polyhedron('V', V_EV_rot), nx);
Xi_true_LV = embedInFullDim(Polyhedron('V', V_LV_rot), nx);

% Xi_true_EV  = Rotation_EV*Polyhedron([-0.06 -0.015;0.06 -0.015; 0.01 0.025; -0.01 0.025]); 
% Xi_true_LV  = Rotation_LV*Polyhedron([-0.06 0.015;0.06 0.015; 0.01 -0.025; -0.01 -0.025]); 

%%% CHANGE! - EMBED INTO HIGHER DIMENSIONS
% Xi_true_EV = embedInFullDim(Xi_true_EV, nx);
% Xi_true_LV = embedInFullDim(Xi_true_LV, nx);

disp(['Dimension of Xi_true_EV: ', num2str(Xi_true_EV.Dim)]);
disp(['Dimension of Xi_true_LV: ', num2str(Xi_true_LV.Dim)]);

min_u_LV    = -1/20;  
max_u_LV    = 1/16;   

%%% CHANGE!
% U_true_LV   = Polyhedron([1/max_u_LV; 1/min_u_LV], [1; 1]);
U_true_LV = Polyhedron([1; -1], [max_u_LV; -min_u_LV]);

W_true      = Xi_true_EV + (-1*Xi_true_LV) + (-B*U_true_LV);  
disp(['Dimension W_true: ', num2str(W_true.Dim)]);


% Above: Construct the true uncertainty set of EV and LV
%%% CHANGE!
W = embedInFullDim(Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2]), nx);

% W           = Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2]);  
% Xi_true_EV  = minHRep(Xi_true_EV); 
% Xi_true_LV  = minHRep(Xi_true_LV); 
W_true      = minHRep(W_true);         
W           = minHRep(W);                   
samples_opt = W_true.V;           
[N_sam_opt, ~] = size(samples_opt);

samples_opt = W_true.V;  % should be n_samples x nx
disp(size(samples_opt)); % should print [num_samples, 3]


distance_error = 15; % constraint on the relative distance between EV and LV
acc_bound      = 2;  % constraint on acceleration of EV

%%% CHANGE! - CONSTRAINTS
F = zeros(nc, nx);
F(1,1) = 1/distance_error;
F(2,1) = -1/distance_error;

G = zeros(nc, nu);
G(3,end) = 1/acc_bound;
G(4,end) = -1/acc_bound;

%G = [0; 0; 1/acc_bound; -1/acc_bound];

% F     = [1/distance_error 0 0; -1/distance_error 0 0; 0 0 0; 0 0 0];       %%% CHANGE! 
% G     = [0; 0; 1/acc_bound; -1/acc_bound];  

%%% CHANGE! - DESIRED STATE
x_des = zeros(nx, 1);
x_des(1) = -35;
% x_des = [-35; 0];    % desired state: relative long. distance and relative long. speed

% nx    = 3;  %%% CHANGE!
% nu    = 1; 
% nc    = 4; 

E     = zeros(1, N);
E(1)  = 1; 
diag_ele = ones(N, 1);
for i    = 1:1:N
    diag_ele(i) = B'*Px*B + R;
end
Pc  = diag(diag_ele); 
M   = diag(ones(N-1, 1), 1); 
Psi = cell(2, 2);
Psi{1, 1} = Phi;
Psi{1, 2} = B*E;
Psi{2, 1} = zeros(nu*N, nx);
Psi{2, 2} = M;
Psi   = cell2mat(Psi);
F_bar = [F + G*K G*E]; 

opts_ini_set.B          = B;
opts_ini_set.N_pre_sam  = 5; 
opts_ini_set.N_sam_opt  = N_sam_opt;
opts_ini_set.W          = W;
opts_ini_set.min_u_LV   = min_u_LV;
opts_ini_set.max_u_LV   = max_u_LV;
opts_ini_set.Xi_true_EV = Xi_true_EV;
opts_ini_set.Xi_true_LV = Xi_true_LV;
opts_ini_set.nx         = nx;
IniSet                  = InitialSetComputation(opts_ini_set);
[alpha_opt, v_opt]      = IniSet.solve_opt(samples_opt');
W_hat_opt               = (1 - alpha_opt)*v_opt + alpha_opt*W; 

%%% CHANGE! convert uncertainty set from ConvexSet (YALMIP) to MPT
%%% Polyhedron
W_hat_opt = Polyhedron('A', W_hat_opt.A, 'b', W_hat_opt.b);

%%% CHANGE! to make it faster
epsilon = 1e-2;

S_hat_opt        = MRPISet(Phi, W_hat_opt, epsilon);
S_true           = MRPISet(Phi, W_true, epsilon); 
S                = MRPISet(Phi, W, epsilon); 
S_hat_opt        = minHRep(S_hat_opt);
S_true           = minHRep(S_true);      
S                = minHRep(S);                
num_half_space_S = length(S.b); 

opts_feasible_region.S         = S;
opts_feasible_region.S_true    = S_true;
opts_feasible_region.S_hat_opt = S_hat_opt;
opts_feasible_region.F         = F;
opts_feasible_region.G         = G;
opts_feasible_region.K         = K;
opts_feasible_region.F_bar     = F_bar;
opts_feasible_region.Psi       = Psi;
opts_feasible_region.nc        = nc;
opts_feasible_region.nx        = nx;
opts_feasible_region.nu        = nu;
opts_feasible_region.N         = N;
opts_feasible_region.Ps        = ones(nx + N*nu, nx + N*nu);
opts_feasible_region.num_half_space_S = num_half_space_S;
IniFeaReg                             = ComputeFeasibleRegion(opts_feasible_region);
[F_N_RMPC, hs_RMPC, Nu_RMPC]          = IniFeaReg.ComFeasibleRegion_RMPC( ); 
[F_N_True, hs_True]       = IniFeaReg.ComFeasibleRegion_True( ); 
[F_N_Hat_Opt, hs_Hat_Opt] = IniFeaReg.ComFeasibleRegion_UQOPT( ); 
opts_feasible_region.hs   = hs_RMPC;

opts_Car.A          = A;
opts_Car.B          = B;
opts_Car.min_u_LV   = min_u_LV;
opts_Car.max_u_LV   = max_u_LV;
opts_Car.Xi_true_EV = Xi_true_EV;
opts_Car.Xi_true_LV = Xi_true_LV;

opts_RMPC.N     = N;
opts_RMPC.nx    = nx;
opts_RMPC.nu    = nu;
opts_RMPC.nc    = nc;
opts_RMPC.Nu    = Nu_RMPC;
opts_RMPC.A     = A;
opts_RMPC.B     = B;
opts_RMPC.F     = F;
opts_RMPC.G     = G;
opts_RMPC.K     = K;
opts_RMPC.Px    = Px;
opts_RMPC.Pc    = Pc;
opts_RMPC.F_bar = F_bar;
opts_RMPC.Psi   = Psi;
opts_RMPC.Phi   = Phi;
opts_RMPC.E     = E;
opts_RMPC.W     = W;
opts_RMPC.S     = S;
opts_RMPC.num_half_space_S = num_half_space_S;
opts_RMPC.hs               = hs_RMPC;

opts_UQMPC    = opts_RMPC;
opts_UQMPC.Nu = Nu_RMPC; % here the UQ-RMPC does not update \nu at every time step

% Computed Ps
%%%% ATTENTION: YOU NEED sdpt3 solver with Yalmip here
H = [ ];
h = [ ];
for i = 0:1:Nu_RMPC
    H = [H; F_bar*(Psi^i)];
    h = [h; 1 - hs_RMPC];
end
Poly        = Polyhedron(H, h);
Poly        = minHRep(Poly);
V           = Poly.V;
[m, ~]      = size(V); 
Ps          = sdpvar(nx + N*nu, nx + N*nu); 
Constraints = [ ];
for i = 1:1:m
    Constraints = [Constraints, norm(Ps*V(i, :)', 2) <= 1]; 
end
Objective               = -logdet(Ps); 
sol                     = optimize(Constraints, Objective);
Ps                      = value(Ps);
opts_UQMPC.Ps           = Ps;
opts_feasible_region.Ps = Ps;

% save the parameters
parameters.opts_ini_set         = opts_ini_set;
parameters.opts_feasible_region = opts_feasible_region;
parameters.opts_Car             = opts_Car;
parameters.opts_RMPC            = opts_RMPC;
parameters.opts_UQMPC           = opts_UQMPC;
parameters.W_hat_opt            = W_hat_opt;
parameters.F_N_RMPC             = F_N_RMPC;
parameters.F_N_RMPC             = F_N_RMPC;
parameters.F_N_True             = F_N_True;
parameters.F_N_Hat_Opt          = F_N_Hat_Opt;
parameters.x_des                = x_des;
parameters.distance_error       = distance_error;
parameters.acc_bound            = acc_bound;
parameters.T                    = T;
parameters.W_true               = W_true;
save('parameters_NEW.mat', 'parameters')

function P_full = embedInFullDim(P, nx)
    % Pads the vertices of a lower-dimensional polyhedron with zeros
    V = P.V;
    V_full = [V, zeros(size(V,1), nx - size(V,2))];
    P_full = Polyhedron('V', V_full);
end

