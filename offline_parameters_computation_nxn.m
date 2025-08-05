% General Parameters
addpath('Functions_general');
clear
close all

main_tic = tic;
stats = struct();
stats.nx = 4; %%% NUMBER OF STATES
stats.nu = 1; %%% NUMBER OF INPUTS
stats.N = 8;
stats.start_time = datetime('now');

nx = stats.nx;
nu = stats.nu;
N = stats.N;

% --- System Setup ---
T = 0.5;

A = eye(nx)*0.5;
for i = 1:nx-1
    A(i, i+1) = T*0.5;
end
B = zeros(nx, nu);
B(end) = T*0.5;

Q = eye(nx)*1.0;
R = 0.01;

[K, Px]  = dlqr(A, B, Q, R);
K        = -K;
Phi      = A + B*K;

% Below: Construct the true uncertainty set of EV and LV
if nx == 2
    theta_EV = pi/40;           
    theta_LV = pi/40;           
    Rotation_EV = [cos(theta_EV) -sin(theta_EV); sin(theta_EV) cos(theta_EV)]; 
    Rotation_LV = [cos(theta_LV) -sin(theta_LV); sin(theta_LV) cos(theta_LV)]; 
    
    Xi_true_EV  = Rotation_EV*Polyhedron([-0.06 -0.015; 0.06 -0.015; 0.01 0.025; -0.01 0.025]); 
    Xi_true_LV  = Rotation_LV*Polyhedron([-0.06 0.015; 0.06 0.015; 0.01 -0.025; -0.01 -0.025]);

elseif nx == 3
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
        [hex_r_LV_x*cos(hex_angle_LV)', hex_r_LV_y*sin(hex_angle_LV)', z_low*ones(6,1)];
        [hex_r_LV_x*cos(hex_angle_LV)', hex_r_LV_y*sin(hex_angle_LV)', z_high*ones(6,1)]
    ];

    V_EV_rot = (Rotation_EV * V_EV_3D')';
    V_LV_rot = (Rotation_LV * V_LV_3D')';

    Xi_true_EV = Polyhedron('V', V_EV_rot);
    Xi_true_LV = Polyhedron('V', V_LV_rot);

elseif nx == 4
    % Define a 4D cuboid for EV and LV uncertainty sets (box, axis-aligned)
    cuboid_half_lengths_EV = [0.015, 0.01, 0.008, 0.012]; % example values
    cuboid_half_lengths_LV = [0.018, 0.009, 0.007, 0.011]; % example values

    % Generate all 16 vertices for the 4D cuboid
    V_EV_4D = [];
    V_LV_4D = [];
    for i = 0:15
        bits = bitget(i, 4:-1:1)*2 - 1; % [-1,1] for each dimension
        V_EV_4D = [V_EV_4D; bits .* cuboid_half_lengths_EV];
        V_LV_4D = [V_LV_4D; bits .* cuboid_half_lengths_LV];
    end

    % Add a small random rotation for realism
    rng(42); % for reproducibility
    angle = pi/18; % small rotation angle
    % Generate a random orthogonal matrix close to identity
    [Q_EV,~] = qr(eye(4) + angle*randn(4));
    [Q_LV,~] = qr(eye(4) + angle*randn(4));
    V_EV_4D_rot = (Q_EV * V_EV_4D')';
    V_LV_4D_rot = (Q_LV * V_LV_4D')';

    Xi_true_EV = Polyhedron('V', V_EV_4D_rot);
    Xi_true_LV = Polyhedron('V', V_LV_4D_rot);

else
    % For nx > 4, use robust box embedding as before
    theta_EV = pi/12;
    theta_LV = pi/12;
    Rotation_EV = [cos(theta_EV) -sin(theta_EV) 0; sin(theta_EV) cos(theta_EV) 0; 0 0 1];
    Rotation_LV = [cos(theta_LV) -sin(theta_LV) 0; sin(theta_LV) cos(theta_LV) 0; 0 0 1];

    hex_angle = (0:5)*pi/3;
    hex_r_x = 0.01;
    hex_r_y = 0.01;
    z_low = -0.01; z_high = 0.01;
    V_EV_3D = [
        [hex_r_x*cos(hex_angle)', hex_r_y*sin(hex_angle)', z_low*ones(6,1)];
        [hex_r_x*cos(hex_angle)', hex_r_y*sin(hex_angle)', z_high*ones(6,1)]
    ];
    hex_angle_LV = (0:5)*pi/3;
    hex_r_LV_x = 0.01;
    hex_r_LV_y = 0.01;
    V_LV_3D = [
        [hex_r_LV_x*cos(hex_angle_LV)', hex_r_LV_y*sin(hex_angle_LV)', z_low*ones(6,1)];
        [hex_r_LV_x*cos(hex_angle_LV)', hex_r_LV_y*sin(hex_angle_LV)', z_high*ones(6,1)]
    ];

    V_EV_rot = (Rotation_EV * V_EV_3D')';
    V_LV_rot = (Rotation_LV * V_LV_3D')';

    Xi_true_EV = embedInFullDimBox(Polyhedron('V', V_EV_rot), nx, 0.05)*0.5;
    Xi_true_LV = embedInFullDimBox(Polyhedron('V', V_LV_rot), nx, 0.05)*0.5;
end

disp(['Dimension of Xi_true_EV: ', num2str(Xi_true_EV.Dim)]);
disp(['Dimension of Xi_true_LV: ', num2str(Xi_true_LV.Dim)]);

min_u_LV    = -1/20;
max_u_LV    = 1/16;
U_true_LV   = Polyhedron([1/max_u_LV; -1/min_u_LV], [1; -1]);

% --- Robust box embedding function ---
function P_full = embedInFullDimBox(P, nx, pad_value)
    % Robustly embeds a lower-dimensional polyhedron P into nx dimensions by creating a box in new dimensions
    % pad_value: size of the box in new dimensions (default 0.01)
    if nargin < 3
        pad_value = 0.01;
    end
    V = P.V;
    d = size(V,2);
    if nx == d
        P_full = Polyhedron('V', V);
        return;
    end
    % For each vertex, create all combinations in new dimensions
    n_new = nx - d;
    box_corners = dec2bin(0:2^n_new-1) - '0'; % binary corners
    box_corners = pad_value * (2*box_corners - 1); % [-pad_value, pad_value] for each new dim
    V_full = [];
    for i = 1:size(V,1)
        v_base = V(i,:);
        v_box = [repmat(v_base, size(box_corners,1), 1), box_corners];
        V_full = [V_full; v_box];
    end
    P_full = Polyhedron('V', V_full);
end

% U_true_LV = Polyhedron([1; -1], [max_u_LV; -min_u_LV]);

W_true      = Xi_true_EV + (-1*Xi_true_LV) + (-B*U_true_LV);
disp(['Dimension W_true: ', num2str(W_true.Dim)]);
disp('First 10 rows of W_true.V:');
disp(W_true.V(1:min(10, size(W_true.V,1)), :));

% After constructing W_true, ensure it is full-dimensional
if rank(W_true.V - W_true.V(1,:)) < nx
    warning('W_true is degenerate! Adding small random noise to vertices.');
    W_true.V = W_true.V + 1e-4*randn(size(W_true.V));
    W_true = Polyhedron('V', W_true.V);
end

% DEFINING W
if nx == 2
    W = Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2]);
elseif nx == 3
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
    ]*1.5; % Scale factor for robustness, adjust as needed
    W_3D = Polyhedron('V', V_cuboid);    
    W = W_3D;
elseif nx == 4
    % 4D cuboid with 16 vertices
    cuboid_half_lengths = [0.12, 0.08, 0.06, 0.09]; % example values, adjust as needed
    V_cuboid_4D = [];
    for i = 0:15
        bits = bitget(i, 4:-1:1)*2 - 1; % [-1,1] for each dimension
        V_cuboid_4D = [V_cuboid_4D; bits .* cuboid_half_lengths];
    end
    W = Polyhedron('V', V_cuboid_4D);
else
    % For nx >= 5, use robust box embedding
    V_cuboid = [
        -0.12, -0.08, -0.06;
         0.12, -0.08, -0.06;
         0.12,  0.08, -0.06;
        -0.12,  0.08, -0.06;
        -0.12, -0.08,  0.06;
         0.12, -0.08,  0.06;
         0.12,  0.08,  0.06;
        -0.12,  0.08,  0.06
    ] * 0.05;
    W_3D = Polyhedron('V', V_cuboid);
    W = embedInFullDimBox(W_3D, nx, 0.01);
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
elseif nx == 3
    nc_state = 2*nx;
    nc_input = 2;
    nc = nc_state + nc_input;
    F = zeros(nc, nx);
    G = zeros(nc, nu);
    % State constraints (axis-aligned)
    distance_error = 15; % distance error in meters
    state_bounds = distance_error * ones(nx,1);
    acc_bound = 4;
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
elseif nx == 4
    % Use simple axis-aligned constraints like in 3D case
    nc_state = 2*nx;  % Just like 3D case: two constraints per dimension
    nc_input = 2;
    nc = nc_state + nc_input;
    
    % Define constraints
    distance_error = 60;  % Keep same scale as before
    acc_bound = 60;
    state_bounds = distance_error * ones(nx,1);
    
    % State constraints (axis-aligned only)
    F = zeros(nc, nx);
    G = zeros(nc, nu);
    
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
else
    % Use previous axis-aligned method for nx > 4
    nc_state = 2*nx;
    nc_input = 2;
    nc = nc_state + nc_input;
    F = zeros(nc, nx);
    G = zeros(nc, nu);
    distance_error = 60;
    state_bounds = distance_error * ones(nx,1);
    acc_bound = 60;
    for i = 1:nx
        F(2*i-1, i) = 1/state_bounds(i);
        F(2*i,   i) = -1/state_bounds(i);
    end
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

% --- PATCH: Ensure F_bar construction is robust and compatible ---
disp(['Size F: ', mat2str(size(F))]);
disp(['Size G: ', mat2str(size(G))]);
disp(['Size K: ', mat2str(size(K))]);
disp(['Size E: ', mat2str(size(E))]);

% F is (nc, nx), G is (nc, nu), K is (nu, nx)
% G*K is (nc, nx)
% E should be (nu, 1) for G*E to be (nc, 1)
if size(E,1) ~= nu
    E = ones(nu, 1); % ensure E is column vector of length nu
end

F_bar = [F + G*K, G*E];

disp(['Size F_bar: ', mat2str(size(F_bar))]);

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
    % Diagnostics: check rank and variance
    disp(['Rank of W_hat_opt vertices: ', num2str(rank(W_hat_opt.V - W_hat_opt.V(1,:)))]);
    disp(['Variance of W_hat_opt vertices in each dimension:']);
    disp(var(W_hat_opt.V,0,1));
    if rank(W_hat_opt.V - W_hat_opt.V(1,:)) < nx
        warning('W_hat_opt is degenerate or nearly so! Adding small random noise to vertices.');
        W_hat_opt.V = W_hat_opt.V + 1e-4*randn(size(W_hat_opt.V));
        W_hat_opt = Polyhedron('V', W_hat_opt.V);
        disp('W_hat_opt Polyhedron after noise:');
        disp(W_hat_opt.V);
        disp(['Rank after noise: ', num2str(rank(W_hat_opt.V - W_hat_opt.V(1,:)))]);
    end
catch ME
    disp('Failed to create W_hat_opt Polyhedron!');
    disp(ME.message);
    error('Aborting due to W_hat_opt Polyhedron creation failure.');
end

epsilon = 0.1;

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
max_vertices = 40;
if size(W_true.V,1) > max_vertices
    warning('Reducing W_true vertices for MRPISet...');
    idx = randperm(size(W_true.V,1), max_vertices);
    W_true = Polyhedron('V', W_true.V(idx,:));
end

t_mrpi_strue = tic;
try
    disp('Calling MRPISet for S_true...');
    S_true = MRPISet(Phi, W_true, epsilon);
    try
        % --- FIX for CDD Error: Rescale before minHRep ---
        if S_true.hasVRep && ~isempty(S_true.V)
            % 1. Find the scaling factor and center
            V = S_true.V;
            center = mean(V, 1);
            max_abs_coord = max(abs(V - center), [], 'all');
            
            if max_abs_coord > 1e-6 % Avoid division by zero
                scale_factor = 1 / max_abs_coord;

                % 2. Scale and center the polyhedron
                S_true_scaled = Polyhedron('V', (V - center) * scale_factor);

                % 3. Compute minimal representation on the scaled version
                S_true_scaled.minHRep();

                % 4. Scale back to the original size and position
                A_scaled = S_true_scaled.A;
                b_scaled = S_true_scaled.b;
                S_true = Polyhedron('A', A_scaled, 'b', b_scaled + A_scaled * center');
            else
                 S_true.minHRep(); % No scaling needed if already centered at zero
            end
        else
            S_true = minHRep(S_true); % Use original call if no V-rep
        end
    catch ME_minHRep
        if contains(ME_minHRep.message, 'CDD')
            warning('minHRep failed for S_true with a CDD error. Approximating S_true with its bounding box.');
            if S_true.hasVRep
                V = S_true.V;
                min_v = min(V, [], 1);
                max_v = max(V, [], 1);
                S_true = Polyhedron('lb', min_v, 'ub', max_v);
            else
                warning('Cannot create bounding box for S_true as V-representation is missing. Re-throwing error.');
                rethrow(ME_minHRep);
            end
        else
            rethrow(ME_minHRep);
        end
    end
    stats.t_mrpi_strue = toc(t_mrpi_strue);
    stats.S_true_num_vertices = size(S_true.V,1);
    stats.S_true_num_constraints = size(S_true.A,1);
    stats.S_true_isBounded = S_true.isBounded;
    disp('S_true vertices:');
    disp('Number of vertices:');
    disp(size(S_true.V,1));
    disp('Max abs vertex value:');
    disp(max(abs(S_true.V),[],'all'));
    % --- Diagnostics for S_true boundedness and rank ---
    %disp('S_true constraint matrix (A):');
    %disp(S_true.A);
    %disp('S_true constraint vector (b):');
    %disp(S_true.b);
    disp(['Rank of S_true.A: ', num2str(rank(S_true.A))]);
    if rank(S_true.A) < nx
        warning('S_true.A does not constrain all directions! S_true will be unbounded.');
    end
    disp(['Rank of S_true vertices: ', num2str(rank(S_true.V - S_true.V(1,:)))]);

    disp('Is S_true bounded?');
    disp(S_true.isBounded);
    % --- PATCH: Check boundedness and rank, fallback if needed ---
    if ~S_true.isBounded
        warning('S_true is not bounded! Fallback to small box polytope.');
        S_true = construct_box_polytope(nx, 0.05, 1e-4);
        S_true = minHRep(S_true);
        stats.S_true_isBounded = S_true.isBounded;
        disp('Fallback S_true vertices:');
        disp(S_true.V);
        disp('Is fallback S_true bounded?');
        disp(S_true.isBounded);
    end
    if rank(S_true.V - S_true.V(1,:)) < nx
        warning('S_true is degenerate or nearly so! Adding small random noise to vertices.');
        S_true.V = S_true.V + 1e-4*randn(size(S_true.V));
        S_true = Polyhedron('V', S_true.V);
        S_true = minHRep(S_true);
        disp('S_true Polyhedron after noise:');
        disp(S_true.V);
        disp(['Rank after noise: ', num2str(rank(S_true.V - S_true.V(1,:)))]);
    end
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

% After MRPISet for S_true, ensure boundedness and full rank
if ~S_true.isBounded
    warning('S_true is not bounded! Fallback to small box polytope.');
    S_true = construct_box_polytope(nx, 0.05, 1e-4);
    S_true = minHRep(S_true);
end
if rank(S_true.V - S_true.V(1,:)) < nx
    warning('S_true is degenerate or nearly so! Adding small random noise to vertices.');
    S_true.V = S_true.V + 1e-4*randn(size(S_true.V));
    S_true = Polyhedron('V', S_true.V);
    S_true = minHRep(S_true);
end

% --- Timing: MRPI S ---
t_mrpi_s = tic;
try
    disp('Calling MRPISet for S...');
    S = MRPISet(Phi, W, epsilon);
    try
        % --- FIX for CDD Error: Rescale before minHRep ---
        if S.hasVRep && ~isempty(S.V)
            % 1. Find the scaling factor and center
            V = S.V;
            center = mean(V, 1);
            max_abs_coord = max(abs(V - center), [], 'all');

            if max_abs_coord > 1e-6 % Avoid division by zero
                scale_factor = 1 / max_abs_coord;

                % 2. Scale and center the polyhedron
                S_scaled = Polyhedron('V', (V - center) * scale_factor);

                % 3. Compute minimal representation on the scaled version
                S_scaled.minHRep();

                % 4. Scale back to the original size and position
                A_scaled = S_scaled.A;
                b_scaled = S_scaled.b;
                S = Polyhedron('A', A_scaled, 'b', b_scaled + A_scaled * center');
            else
                S.minHRep(); % No scaling needed if already centered at zero
            end
        else
            S = minHRep(S); % Use original call if no V-rep
        end
    catch ME_minHRep
        if contains(ME_minHRep.message, 'CDD')
            warning('minHRep failed for S with a CDD error. Approximating S with its bounding box.');
            if S.hasVRep
                V = S.V;
                min_v = min(V, [], 1);
                max_v = max(V, [], 1);
                S = Polyhedron('lb', min_v, 'ub', max_v);
            else
                warning('Cannot create bounding box for S as V-representation is missing. Re-throwing error.');
                rethrow(ME_minHRep);
            end
        else
            rethrow(ME_minHRep);
        end
    end
    stats.t_mrpi_s = toc(t_mrpi_s);
    stats.S_num_vertices = size(S.V,1);
    stats.S_num_constraints = size(S.A,1);
    stats.S_isBounded = S.isBounded;
    disp('S vertices:');
    disp(size(S.V,1));
    disp('Max abs vertex value:');
    disp(max(abs(S.V),[],'all'));
    disp('Is S bounded?');
    disp(S.isBounded);
    % Diagnostic: print W and Phi
    disp('Diagnostic: W vertices:');
    disp(W.V);
    disp('Diagnostic: W max abs vertex value:');
    disp(max(abs(W.V),[],'all'));
    disp('Diagnostic: Phi matrix:');
    disp(Phi);
    disp('Diagnostic: Phi eigenvalues:');
    disp(eig(Phi));
    disp('Diagnostic: epsilon:');
    disp(epsilon);
    % Additional diagnostics for S constraint matrices
    if isprop(S, 'A') && isprop(S, 'b')
        disp('Diagnostic: S.A size:');
        disp(size(S.A));
        disp('Diagnostic: S.b size:');
        disp(size(S.b));
        if isempty(S.A) || isempty(S.b)
            warning('S.A or S.b is empty! MRPI set S may be degenerate or empty.');
        end
    end
catch ME
    stats.t_mrpi_s = toc(t_mrpi_s);
    disp('MRPISet or minHRep failed for S, even after fallback!');
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

% --- Automated post-processing and diagnostics for hs ---
hs_threshold = 1e6;
large_hs_idx = find(abs(hs_RMPC) > hs_threshold);
if ~isempty(large_hs_idx)
    warning('Some hs values are extremely large! Indices:');
    disp(large_hs_idx);
    disp('Large hs values:');
    disp(hs_RMPC(large_hs_idx));
    % Cap large hs values for downstream use (optional)
    hs_RMPC_capped = hs_RMPC;
    hs_RMPC_capped(large_hs_idx) = sign(hs_RMPC(large_hs_idx)) * hs_threshold;
    disp('Capped hs values:');
    disp(hs_RMPC_capped(large_hs_idx));
    % Print condition number diagnostics for F_Com
    if exist('F_N_RMPC','var')
        cond_F = cond(F_N_RMPC);
        disp(['Condition number of F_N_RMPC: ', num2str(cond_F)]);
        if cond_F > 1e8
            warning('F_N_RMPC is ill-conditioned!');
        end
    end
    % Optionally, flag for further relaxation or fallback
    % For full automation, you could retry with inflated W or increased epsilon here
    % For now, just cap and continue
    opts_feasible_region.hs = hs_RMPC_capped;
else
    opts_feasible_region.hs = hs_RMPC;
end

function vol = bounding_box_volume(P)
    % Computes the volume of the axis-aligned bounding box of polyhedron P
    v_min = min(P.V, [], 1);
    v_max = max(P.V, [], 1);
    vol = prod(v_max - v_min);
end

stats.F_N_RMPC_bbox_vol = bounding_box_volume(F_N_RMPC);
stats.F_N_True_bbox_vol = bounding_box_volume(F_N_True);
stats.F_N_Hat_Opt_bbox_vol = bounding_box_volume(F_N_Hat_Opt);

fprintf('F_N_RMPC bounding box volume: %.4g\n', stats.F_N_RMPC_bbox_vol);
fprintf('F_N_True bounding box volume: %.4g\n', stats.F_N_True_bbox_vol);
fprintf('F_N_Hat_Opt bounding box volume: %.4g\n', stats.F_N_Hat_Opt_bbox_vol);


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
save(['stats_FINAL_nx', num2str(nx), feas_str_stats, '_2.mat'], 'stats');

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
save(['parameters_FINAL_nx', num2str(nx), feas_str, '_3.mat'], 'parameters');

fprintf('\nComputation time: %.2f seconds\n', stats.elapsedTime)

% --- Helper Functions ---
%{
function P_full = embedInFullDim(P, nx)
    % Pads the vertices of a lower-dimensional polyhedron with zeros
    V = P.V;
    V_full = [V, zeros(size(V,1), nx - size(V,2))];
    P_full = Polyhedron('V', V_full);
end
%}

function P_full = embedInFullDim(P, nx, pad_value)
    % Pads the vertices of a lower-dimensional polyhedron with pad_value (default 0)
    if nargin < 3
        pad_value = 0;
    end
    V = P.V;
    % If pad_value is a scalar, add small random noise to avoid degeneracy
    if pad_value ~= 0
        noise = 1e-5 * randn(size(V,1), nx - size(V,2));
        V_full = [V, pad_value * ones(size(V,1), nx - size(V,2)) + noise];
    else
        V_full = [V, zeros(size(V,1), nx - size(V,2))];
    end
    if any(~isfinite(V_full(:)))
        error('Embedded vertices contain Inf or NaN!');
    end
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

function Ps = compute_Ps(F_bar, Psi, Nu_RMPC, hs_RMPC, nx, N, nu)
    % Compute the constraint polyhedron for the scenario-based MPC
    disp('compute_Ps:');
    H = [];
    h = [];
    for i = 0:Nu_RMPC
        H = [H; F_bar*(Psi^i)];
        h = [h; 1 - hs_RMPC];
    end
    % Use H-representation directly for constraints
    [m, n] = size(H);
    Ps = sdpvar(nx + N*nu, nx + N*nu);
    Constraints = [];
    for i = 1:m
        % Use the rows of H as directions, or sample feasible points instead
        % Example: norm(Ps*H(i,:)', 2) <= 1; (or sample feasible z)
        Constraints = [Constraints, norm(Ps*H(i,:)', 2) <= 1];
    end
    Objective = -logdet(Ps);
    optimize(Constraints, Objective);
    Ps = value(Ps);
end