% General Parameters
addpath('Functions_general');
clear
close all

main_tic = tic;
stats = struct();
stats.nx = 2; %%% NUMBER OF STATES
stats.nu = 1; %%% NUMBER OF INPUTS
stats.N = 8;
stats.start_time = datetime('now');

nx = stats.nx;
nu = stats.nu;
N = stats.N;

% --- System Setup ---
T = 0.5;

A = eye(nx)*1.0;
for i = 1:nx-1
    A(i, i+1) = T*1.0;
end
B = zeros(nx, nu);
B(end) = T;

Q = eye(nx)*1;
R = 0.1;

[K, Px]  = dlqr(A, B, Q, R);
K        = -K;
Phi      = A + B*K;

% Below: Construct the true uncertainty set of EV and LV
if nx == 2
    theta_EV = pi/40;           
    theta_LV = pi/40;           
    Rotation_EV = [cos(theta_EV) -sin(theta_EV); sin(theta_EV) cos(theta_EV)]; 
    Rotation_LV = [cos(theta_LV) -sin(theta_LV); sin(theta_LV) cos(theta_LV)]; 
    
    Xi_true_EV  = Rotation_EV*Polyhedron([-0.06 -0.015;0.06 -0.015; 0.01 0.025; -0.01 0.025]); 
    Xi_true_LV  = Rotation_LV*Polyhedron([-0.06 0.015;0.06 0.015; 0.01 -0.025; -0.01 -0.025]);

else % 3D+ case
    theta_EV = pi/12;
    theta_LV = pi/12;
    % Define 3D rotation matrix (e.g., around z-axis)
    Rotation_EV = [cos(theta_EV) -sin(theta_EV) 0; sin(theta_EV) cos(theta_EV) 0; 0 0 1];
    Rotation_LV = [cos(theta_LV) -sin(theta_LV) 0; sin(theta_LV) cos(theta_LV) 0; 0 0 1];

    % Define 3D vertices for EV uncertainty (asymmetric hexagon in x-y, extruded in z)
    hex_angle = (0:5)*pi/3; % 6 points for hexagon
    hex_r_x = 0.018; % x-radius (longer)
    hex_r_y = 0.012; % y-radius (shorter)
    z_low = -0.005; z_high = 0.005;
    V_EV_3D = [
        [hex_r_x*cos(hex_angle)', hex_r_y*sin(hex_angle)', z_low*ones(6,1)];
        [hex_r_x*cos(hex_angle)', hex_r_y*sin(hex_angle)', z_high*ones(6,1)]
    ];
    % Define 3D vertices for LV uncertainty (asymmetric hexagon in x-y, extruded in z)
    hex_angle_LV = (0:5)*pi/3;
    hex_r_LV_x = 0.022; % x-radius (longer)
    hex_r_LV_y = 0.010; % y-radius (shorter)
    V_LV_3D = [
        [hex_r_LV_x*cos(hex_angle_LV)', hex_r_LV_y*sin(hex_angle_LV)', z_low*ones(6,1)]*2;
        [hex_r_LV_x*cos(hex_angle_LV)', hex_r_LV_y*sin(hex_angle_LV)', z_high*ones(6,1)]*2
    ];

    V_EV_rot = (Rotation_EV * V_EV_3D')';
    V_LV_rot = (Rotation_LV * V_LV_3D')';

    % Embed the 3D polyhedron into the full nx-dimensional space
    Xi_true_EV = embedInFullDim(Polyhedron('V', V_EV_rot), nx);
    Xi_true_LV = embedInFullDim(Polyhedron('V', V_LV_rot), nx);
end

disp(['Dimension of Xi_true_EV: ', num2str(Xi_true_EV.Dim)]);
disp(['Dimension of Xi_true_LV: ', num2str(Xi_true_LV.Dim)]);

min_u_LV    = -1/20;
max_u_LV    = 1/16;
U_true_LV   = Polyhedron([1/max_u_LV; 1/min_u_LV], [1; 1]);

% U_true_LV = Polyhedron([1; -1], [max_u_LV; -min_u_LV]);

W_true      = Xi_true_EV + (-1*Xi_true_LV) + (-B*U_true_LV);
disp(['Dimension W_true: ', num2str(W_true.Dim)]);

print_set_info('Xi_true_EV', Xi_true_EV);
print_set_info('Xi_true_LV', Xi_true_LV);
print_set_info('W_true', W_true);

% DEFINING W
if nx == 2
    W = Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2]);
else
    % --- Original/previous implementation (commented out) ---
    % W = embedInFullDim(Polyhedron([-0.01 -0.01 0; 0.01 -0.01 0; 0.01 0.01 0; -0.01 0.01 0]), nx);
    % W = embedInFullDim(Polyhedron([-0.4 -0.2 -0.1; 0.01 -0.01 0; 0.01 0.01 0; -0.01 0.01 0]), nx);

    % --- New 3D cuboid definition ---
    % Define a 3D cuboid with 8 vertices
    V_cuboid = [
        -0.12, -0.08, -0.06;
         0.12, -0.08, -0.06;
         0.12,  0.08, -0.06;
        -0.12,  0.08, -0.06;
        -0.12, -0.08,  0.06;
         0.12, -0.08,  0.06;
         0.12,  0.08,  0.06;
        -0.12,  0.08,  0.06
    ];
    W_3D = Polyhedron('V', V_cuboid);
    W = embedInFullDim(W_3D, nx);

    %%%%% fallback
    %W = embedInFullDim(Polyhedron([-0.01 -0.01 0; 0.01 -0.01 0; 0.01 0.01 0; -0.01 0.01 0]), nx);
end

W_true      = minHRep(W_true);         
W           = minHRep(W);                  

samples_opt = W_true.V;
[N_sam_opt, ~] = size(samples_opt);
print_set_info('samples_opt', Polyhedron('V', samples_opt));

if nx == 2
    nc = 4;
    distance_error = 15;
    acc_bound      = 2;
    F = zeros(nc, nx);
    F(1,1) = 1/distance_error;
    F(2,1) = -1/distance_error;
    G = zeros(nc, nu);
    G(3,end) = 1/acc_bound;
    G(4,end) = -1/acc_bound;
else
    nc_state = 2*nx;
    nc_input = 2;
    nc = nc_state + nc_input;
    F = zeros(nc, nx);
    G = zeros(nc, nu);
    % State constraints
    state_bounds = 50 * ones(nx,1);
    acc_bound = 60;
    for i = 1:nx
        F(2*i-1, i) = 1/state_bounds(i);
        F(2*i,   i) = -1/state_bounds(i);
    end
    % Input constraints
    G(2*nx+1,1) = 1/acc_bound;
    G(2*nx+2,1) = -1/acc_bound;
    disp('F matrix:');
    disp(F);
    disp('G matrix:');
    disp(G);
end

x_des = zeros(nx, 1);
x_des(1) = -35;
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

% --- Timing: InitialSetComputation ---
t_ini = tic;
IniSet = InitialSetComputation(opts_ini_set);
[alpha_opt, v_opt] = IniSet.solve_opt(samples_opt');
stats.t_initialset = toc(t_ini);
stats.alpha_opt = alpha_opt;
stats.v_opt = v_opt;
W_hat_opt = (1 - alpha_opt)*v_opt + alpha_opt*W;

try
    W_hat_opt = Polyhedron('A', W_hat_opt.A, 'b', W_hat_opt.b);
    disp('W_hat_opt Polyhedron created. Vertices:');
    disp(W_hat_opt.V);
catch ME
    disp('Failed to create W_hat_opt Polyhedron!');
    disp(ME.message);
    error('Aborting due to W_hat_opt Polyhedron creation failure.');
end

epsilon = 1e-2;

eig_Phi = eig(Phi);
disp('Closed-loop eigenvalues:');
disp(eig_Phi);
if any(abs(eig_Phi) >= 1)
    warning('Closed-loop system is not stable! MRPI set will be unbounded.');
end

% --- Timing: MRPI S_hat_opt ---
t_mrpi_shat = tic;
try
    disp('Calling MRPISet for S_hat_opt...');
    S_hat_opt = MRPISet(Phi, W_hat_opt, epsilon);
    S_hat_opt = minHRep(S_hat_opt);
    stats.t_mrpi_shat = toc(t_mrpi_shat);
    stats.S_hat_opt_num_vertices = size(S_hat_opt.V,1);
    stats.S_hat_opt_num_constraints = size(S_hat_opt.A,1);
    stats.S_hat_opt_isBounded = S_hat_opt.isBounded;
    disp('S_hat_opt vertices:');
    disp(S_hat_opt.V);
    disp('Number of vertices:');
    disp(size(S_hat_opt.V,1));
    disp('Max abs vertex value:');
    disp(max(abs(S_hat_opt.V),[],'all'));
    disp('Is S_hat_opt bounded?');
    disp(S_hat_opt.isBounded);
catch ME
    stats.t_mrpi_shat = toc(t_mrpi_shat);
    disp('MRPISet or minHRep failed for S_hat_opt!');
    disp(ME.message);
    disp('Diagnostics:');
    disp('Phi:');
    disp(Phi);
    disp('W_hat_opt:');
    disp(W_hat_opt);
    error('Aborting due to unbounded or degenerate MRPI set S_hat_opt.');
end

% --- Timing: MRPI S_true ---
t_mrpi_strue = tic;
try
    disp('Calling MRPISet for S_true...');
    S_true = MRPISet(Phi, W_true, epsilon);
    S_true = minHRep(S_true);
    stats.t_mrpi_strue = toc(t_mrpi_strue);
    stats.S_true_num_vertices = size(S_true.V,1);
    stats.S_true_num_constraints = size(S_true.A,1);
    stats.S_true_isBounded = S_true.isBounded;
    disp('S_true vertices:');
    disp(S_true.V);
catch ME
    stats.t_mrpi_strue = toc(t_mrpi_strue);
    disp('MRPISet or minHRep failed for S_true!');
    disp(ME.message);
    disp('Diagnostics:');
    disp('Phi:');
    disp(Phi);
    disp('W_true:');
    disp(W_true);
    error('Aborting due to unbounded or degenerate S_true.');
end

% --- Timing: MRPI S ---
t_mrpi_s = tic;
try
    disp('Calling MRPISet for S...');
    S = MRPISet(Phi, W, epsilon);
    S = minHRep(S);
    stats.t_mrpi_s = toc(t_mrpi_s);
    stats.S_num_vertices = size(S.V,1);
    stats.S_num_constraints = size(S.A,1);
    stats.S_isBounded = S.isBounded;
    disp('S vertices:');
    disp(S.V);
catch ME
    stats.t_mrpi_s = toc(t_mrpi_s);
    disp('MRPISet or minHRep failed for S!');
    disp(ME.message);
    disp('Diagnostics:');
    disp('Phi:');
    disp(Phi);
    disp('W:');
    disp(W);
    error('Aborting due to unbounded or degenerate S.');
end

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

t_feareg = tic;
IniFeaReg = ComputeFeasibleRegion(opts_feasible_region);
[F_N_RMPC, hs_RMPC, Nu_RMPC] = IniFeaReg.ComFeasibleRegion_RMPC();
[F_N_True, hs_True] = IniFeaReg.ComFeasibleRegion_True();
[F_N_Hat_Opt, hs_Hat_Opt] = IniFeaReg.ComFeasibleRegion_UQOPT();
stats.t_feareg = toc(t_feareg);
stats.F_N_RMPC_size = size(F_N_RMPC);
stats.F_N_True_size = size(F_N_True);
stats.F_N_Hat_Opt_size = size(F_N_Hat_Opt);
opts_feasible_region.hs = hs_RMPC;

opts_Car.A          = A;
opts_Car.B          = B;
opts_Car.min_u_LV   = min_u_LV;
opts_Car.max_u_LV   = max_u_LV;
opts_Car.Xi_true_EV = Xi_true_EV;
opts_Car.Xi_true_LV = Xi_true_LV;

opts_RMPC = struct();
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
opts_UQMPC.Nu = Nu_RMPC;

% --- Timing: SDP ---
t_sdp = tic;
Ps = compute_Ps(F_bar, Psi, Nu_RMPC, hs_RMPC, nx, N, nu);
stats.t_sdp = toc(t_sdp);
opts_UQMPC.Ps = Ps;
opts_feasible_region.Ps = Ps;

% --- Save stats and parameters ---
stats.elapsedTime = toc(main_tic);
feas_str_stats = '';
if exist('Nu_RMPC','var')
    feas_str_stats = ['_feas', num2str(Nu_RMPC/opts_RMPC.N, '%.3f')];
end
save(['stats_FINAL_nx', num2str(nx), feas_str_stats, '.mat'], 'stats');

parameters.opts_ini_set         = opts_ini_set;
parameters.opts_feasible_region = opts_feasible_region;
parameters.opts_Car             = opts_Car;
parameters.opts_RMPC            = opts_RMPC;
parameters.opts_UQMPC           = opts_UQMPC;
parameters.W_hat_opt            = W_hat_opt;
parameters.F_N_RMPC             = F_N_RMPC;
parameters.F_N_True             = F_N_True;
parameters.F_N_Hat_Opt          = F_N_Hat_Opt;
parameters.x_des                = x_des;
parameters.distance_error       = distance_error;
parameters.acc_bound            = acc_bound;
parameters.T                    = T;
parameters.W_true               = W_true;

feas_str = '';
if exist('Nu_RMPC','var')
    feas_str = ['_feas', num2str(Nu_RMPC/opts_RMPC.N, '%.3f')];
end
save(['parameters_FINAL_nx', num2str(nx), feas_str, '.mat'], 'parameters');

fprintf('\nComputation time: %.2f seconds\n', stats.elapsedTime)

% --- Helper Functions ---
function P_full = embedInFullDim(P, nx)
    % Pads the vertices of a lower-dimensional polyhedron with zeros
    V = P.V;
    V_full = [V, zeros(size(V,1), nx - size(V,2))];
    P_full = Polyhedron('V', V_full);
end

function print_set_info(name, P)
    disp(['--- DIAGNOSTIC: ', name, ' ---']);
    disp([name, ' vertices:']); disp(P.V);
    disp([name, ' rank: ', num2str(rank(P.V - P.V(1,:)))]);
    disp(['Variance of ', name, ' vertices in each dimension:']); disp(var(P.V,0,1));
    disp(['--- END DIAGNOSTIC ---']);
end

function W = construct_box_polytope(nx, scaling, noise)
    w_box = 1e-3 * scaling * ones(1, nx);
    n_vertices = 2^nx;
    V = zeros(n_vertices, nx);
    for v = 0:(n_vertices-1)
        bits = bitget(v, nx:-1:1);
        V(v+1, :) = (-w_box) + 2*w_box.*bits;
    end
    V = V + noise*randn(size(V));
    W = Polyhedron('V', V);
end

function [W, success] = try_construct_W(nx, scaling_factors, noise)
    success = false;
    for attempt = 1:length(scaling_factors)
        try
            if nx == 2
                W = Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2] * scaling_factors(attempt));
            else
                W = construct_box_polytope(nx, scaling_factors(attempt), noise);
            end
            W = minHRep(W);
            if rank(W.V - W.V(1,:)) < nx
                warning('W is degenerate or nearly so!');
                continue;
            end
            success = true;
            break;
        catch
            continue;
        end
    end
    if ~success
        error('Failed to create a valid, non-degenerate W after multiple attempts.');
    end
end

function [S, success, W_out] = try_construct_S(Phi, W, scaling_factors, epsilon, nx)
    success = false;
    W_out = W;
    for attempt_s = 1:length(scaling_factors)
        try
            W_scaled = W_out * scaling_factors(attempt_s);
            S_tmp = MRPISet(Phi, W_scaled, epsilon);
            S_tmp = minHRep(S_tmp);
            if rank(S_tmp.V - S_tmp.V(1,:)) < nx
                warning('S is degenerate or nearly so!');
                continue;
            end
            S = S_tmp;
            W_out = W_scaled;
            success = true;
            break;
        catch
            continue;
        end
    end
    if ~success
        error('Failed to create a valid S after multiple scaling attempts.');
    end
end

function Ps = compute_Ps(F_bar, Psi, Nu_RMPC, hs_RMPC, nx, N, nu)
    % Compute the constraint polyhedron for the scenario-based MPC
    H = [];
    h = [];
    for i = 0:Nu_RMPC
        H = [H; F_bar*(Psi^i)];
        h = [h; 1 - hs_RMPC];
    end
    Poly = Polyhedron(H, h);
    Poly = minHRep(Poly);
    V = Poly.V;
    [m, ~] = size(V);
    Ps = sdpvar(nx + N*nu, nx + N*nu);
    Constraints = [];
    for i = 1:m
        Constraints = [Constraints, norm(Ps*V(i, :)', 2) <= 1];
    end
    Objective = -logdet(Ps);
    optimize(Constraints, Objective);
    Ps = value(Ps);
end