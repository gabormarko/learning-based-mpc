% hs_debug_experiment.m
% Standalone script to debug and experiment with hs (robust tube) calculation

addpath('Functions_general');
clear; %clc; 
% close all;

% ==== USER PARAMETERS ====
nx = 3;    % Number of states
nu = 1;    % Number of inputs

% System matrices (edit as needed)
A = eye(nx);
T = 0.5;
for i = 1:nx-1
    A(i, i+1) = T;
end
B = zeros(nx, nu);
B(end) = T;

% LQR gain (edit Q, R as needed)
Q = eye(nx);
R = 0.1;
[K, Px] = dlqr(A, B, Q, R);
K = -K;
Phi = A + B*K;

if nx == 2      %%%% 2D
    nc = 4; % number of constraints
    distance_error = 15; % constraint on the relative distance between EV and LV
    acc_bound      = 2;  % constraint on acceleration of EV

    %%% CHANGE! - CONSTRAINTS
    F = zeros(nc, nx);
    F(1,1) = 1/distance_error;
    F(2,1) = -1/distance_error;

    G = zeros(nc, nu);
    G(3,end) = 1/acc_bound;
    G(4,end) = -1/acc_bound;

    nc_state = 2; % 2 state constraints
else    %%%% 3D+
    % Add 6 state constraints and 2 input constraints
    nc_state = 6; % 2 per state
    nc_input = 2; % 2 for input
    nc = nc_state + nc_input;

    F = zeros(nc, nx);
    G = zeros(nc, nu);

    distance_error = 30;
    velocity_error = 30;
    accel_error    = 30;
    acc_bound      = 10;

    % State constraints
    F(1,1) = 1/distance_error;   % x1 upper
    F(2,1) = -1/distance_error;  % x1 lower
    F(3,2) = 1/velocity_error;   % x2 upper
    F(4,2) = -1/velocity_error;  % x2 lower
    F(5,3) = 1/accel_error;      % x3 upper
    F(6,3) = -1/accel_error;     % x3 lower

    % Input constraints
    G(7,1) = 1/acc_bound;        % u upper
    G(8,1) = -1/acc_bound;       % u lower
end

disp('F matrix:'); disp(F);
disp('G matrix:'); disp(G);

disp('System matrix A:'); disp(A);
disp('System matrix B:'); disp(B);
disp('LQR gain K:'); disp(K);
disp('Closed-loop Phi:'); disp(Phi);

% Uncertainty set (edit as needed)
if nx == 2
    W = Polyhedron([-0.5 -0.2; 0.5 -0.2; 0.5 0.2; -0.5 0.2]);
else
    W = embedInFullDim(Polyhedron([-0.01 -0.01 0; 0.01 -0.01 0; 0.01 0.01 0; -0.01 0.01 0]), nx);
end
disp('Uncertainty set W vertices:'); disp(W.V);

epsilon = 1e-2;
try
    disp('Computing MRPI set S...');
    S = MRPISet(Phi, W, epsilon);
    S = minHRep(S);
    % disp('MRPI set S (A, b):');
    % disp(S.A);
    % disp(S.b);
    num_half_space_S = length(S.b);

    % Only use state constraints for hs calculation
    F_state = F(1:nc_state, :);
    G_state = G(1:nc_state, :);

    % Compute hs (robust tube) for state constraints
    hs = zeros(nc_state, 1);
    disp('Computing hs for each state constraint:');
    for i = 1:nc_state
        % For each constraint, maximize over S: F(i,:)*x <= 1 - hs(i)
        % hs(i) = max_{s in S} F(i,:)*s
        hs(i) = max(F_state(i,:) * S.V');
        %fprintf('  hs(%d) = max(F(%d,:)*S.V) = %g\n', i, i, hs(i));
        % Print the values for each vertex
        %fprintf('    F(%d,:)*S.V = ', i);
        %disp(F_state(i,:) * S.V');
    end

    disp('Computed hs (robust tube):');
    disp(hs);

    % Check feasibility: hs < 1 for all constraints
    if all(hs < 1)
        disp('SUCCESS: hs calculation converged and is feasible for all constraints.');
    else
        disp('WARNING: hs calculation did NOT converge for some constraints (hs >= 1).');
        disp('Try tightening constraints or reducing uncertainty set.');
    end

    % --- Plotting section ---
    figure;
    subplot(1,2,1);
    hold on;
    if nx == 2
        % Plot MRPI set S in 2D
        plot(S, 'color', 'b', 'alpha', 0.2);
        title('MRPI Set S (2D)');
        xlabel('x_1'); ylabel('x_2');
    elseif nx == 3
        % Plot MRPI set S in 3D
        plot(S, 'color', 'b', 'alpha', 0.2);
        title('MRPI Set S (3D)');
        xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
        view(3);
    else
        title('MRPI Set S (plotting only for nx=2 or 3)');
    end
    hold off;

    subplot(1,2,2);
    bar(hs);
    title('hs (robust tube) for each state constraint');
    xlabel('Constraint index'); ylabel('hs value');
    ylim([0, max(1, max(hs)+0.1)]);
    grid on;
    % --- End plotting ---

catch ME
    disp('ERROR during hs calculation:');
    disp(ME.message);
    disp('Stack trace:');
    disp(ME.getReport());
end

% --- LQR gain sweep for aggressive control ---
Q_candidates = {eye(nx), 10*eye(nx), diag([100 10 1]), diag([1000 100 10])};
R_candidates = [0.1, 0.01, 0.001];

best_hs = Inf;
best_Q = Q;
best_R = R;
best_K = K;
best_Phi = Phi;
best_S = [];
best_hs_vec = [];

for qidx = 1:length(Q_candidates)
    for ridx = 1:length(R_candidates)
        Q_try = Q_candidates{qidx};
        R_try = R_candidates(ridx);
        [K_try, ~] = dlqr(A, B, Q_try, R_try);
        K_try = -K_try;
        Phi_try = A + B*K_try;
        S_try = MRPISet(Phi_try, W, epsilon);
        S_try = minHRep(S_try);
        if isempty(S_try.V)
            continue;
        end
        hs_try = zeros(nc_state, 1);
        for i = 1:nc_state
            hs_try(i) = max(F(i,:) * S_try.V');
        end
        max_hs = max(hs_try);
        if max_hs < best_hs
            best_hs = max_hs;
            best_Q = Q_try;
            best_R = R_try;
            best_K = K_try;
            best_Phi = Phi_try;
            best_S = S_try;
            best_hs_vec = hs_try;
        end
    end
end

fprintf('\nBest LQR tuning found: max(hs) = %g\n', best_hs);
disp('Best Q:'); disp(best_Q);
disp('Best R:'); disp(best_R);
disp('Best K:'); disp(best_K);
disp('Best Phi:'); disp(best_Phi);

if best_hs < 1
    disp('SUCCESS: With aggressive LQR tuning, robust tube fits inside constraints.');
else
    disp('WARNING: Even with aggressive LQR tuning, robust tube does not fit. Try loosening constraints or reducing W further.');
end

% Optionally, plot the best MRPI set and hs values
figure;
subplot(1,2,1);
hold on;
if nx == 2
    plot(best_S, 'color', 'g', 'alpha', 0.2);
    title('Best MRPI Set S (2D, LQR sweep)');
    xlabel('x_1'); ylabel('x_2');
elseif nx == 3
    plot(best_S, 'color', 'g', 'alpha', 0.2);
    title('Best MRPI Set S (3D, LQR sweep)');
    xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
    view(3);
else
    title('Best MRPI Set S (plotting only for nx=2 or 3)');
end
hold off;

subplot(1,2,2);
bar(best_hs_vec);
title('Best hs (robust tube) for each state constraint');
xlabel('Constraint index'); ylabel('hs value');
ylim([0, max(1, max(best_hs_vec)+0.1)]);
grid on;
% --- End LQR sweep plotting ---

% --- End of script ---

function P_full = embedInFullDim(P, nx)
    % Pads the vertices of a lower-dimensional polyhedron with zeros
    V = P.V;
    V_full = [V, zeros(size(V,1), nx - size(V,2))];
    P_full = Polyhedron('V', V_full);
end