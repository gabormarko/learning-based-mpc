% Comparison between UQ-RMPC and SCMPC
clc
clear
close all
%load('parameters_NEW.mat');
%load('parameters_2x2.mat');
load('parameters_3x3.mat');
%load('parameters_4x4.mat');

N_Samples = 400;
% Assign parameters for UQ-RMPC
opts_ini_set = parameters.opts_ini_set;
opts_Car = parameters.opts_Car;
opts_UQMPC = parameters.opts_UQMPC;
F_N_True = parameters.F_N_True;
x_des = parameters.x_des;
nx = opts_UQMPC.nx;
nu = opts_UQMPC.nu;
nc = opts_UQMPC.nc;
A = opts_UQMPC.A;
B = opts_UQMPC.B;
T = parameters.T;
% Assign parameters for SCMPC
opts_SCMPC =  parameters.opts_Car;
opts_SCMPC.N_SC = N_Samples/opts_UQMPC.N;
% Generalize Q and R for any nx, nu
opts_SCMPC.Q = eye(nx);
if nu == 1
    opts_SCMPC.R = 0.1;
else
    opts_SCMPC.R = eye(nu)*0.1;
end
opts_SCMPC.N = opts_UQMPC.N;
opts_SCMPC.F = opts_UQMPC.F;
opts_SCMPC.G = opts_UQMPC.G;
opts_SCMPC.nu = opts_UQMPC.nu;
opts_SCMPC.nx = opts_UQMPC.nx;
opts_SCMPC.Xi_true_EV_A = opts_SCMPC.Xi_true_EV.A;
opts_SCMPC.Xi_true_EV_b = opts_SCMPC.Xi_true_EV.b;
opts_SCMPC.Xi_true_LV_A = opts_SCMPC.Xi_true_LV.A;
opts_SCMPC.Xi_true_LV_b = opts_SCMPC.Xi_true_LV.b;

% call MPCs and Car modeling
UQRobustMPC = UQRobustMPC(opts_UQMPC); 
SCMPC       = ScenarioMPC(opts_SCMPC);
Car         = ModelingCar(opts_Car); 

% change the number of samples for initial quantification
opts_ini_set.N_pre_sam = N_Samples;
IniSet = InitialSetComputation(opts_ini_set);
[alpha_ini, v_ini, samples] = IniSet.solve();

% Online iteraction
K_N = 50;          % the MPC iterates 100 steps online

% Generalize initial states for any nx
x_LV_0 = zeros(nx,1);
x_LV_0(1:min(2,nx)) = [100; 10];
x_RM_0 = zeros(nx,1);
x_RM_0(1:min(2,nx)) = [-12; 5];
% Optionally, allow user to set other initial values for nx>2
x_EV_0 = x_LV_0 + x_des + x_RM_0; % initial state of EV

State_LV = ones(nx, K_N + 1); % state of the LV
State_LV(:, 1) = x_LV_0;   
Input_LV = ones(nx, K_N);     % input of the LV, the input is control + disturbance

Noise_EV = ones(nx, K_N);           % disturbance of the EV
Control_EV_UQMPC = ones(nu, K_N);   % control input of EV
State_EV_UQMPC = ones(nx, K_N + 1); % state of the EV
State_RM_UQMPC = ones(nx, K_N + 1); % relative state
Nominal_RM_State_UQMPC = ones(nx, K_N + 1); % nominal relative state
Hs_UQMPC       = ones(nc, K_N);             % value of hs vector
State_EV_UQMPC(:, 1) = x_EV_0;
State_RM_UQMPC(:, 1) = x_EV_0 - x_LV_0 - x_des;
Nominal_RM_State_UQMPC(:, 1) = x_RM_0;


Alpha = ones(1, K_N + 1);
S_Hat = cell(K_N, 1);
W_Hat = cell(K_N, 1);
Alpha(:, 1) = alpha_ini;
% V can be a vector or Polyhedron, so use cell array for generality
V = cell(1, K_N + 1);
V{1} = v_ini;

Control_EV_SCMPC = ones(nu, K_N);   % control input of EV
State_EV_SCMPC = ones(nx, K_N + 1); % state of the EV
State_RM_SCMPC = ones(nx, K_N + 1); % relative state
State_EV_SCMPC(:, 1) = x_EV_0;
State_RM_SCMPC(:, 1) = x_EV_0 - x_LV_0 - x_des;

Time_UQMPC = ones(K_N, 1);
Time_SCMPC = ones(K_N, 1);

for k = 1:K_N
    [xi_LV_k, x_LV_k_next] = Car.LVModeling(State_LV(:, k));
    State_LV(:, k + 1) = x_LV_k_next;
    Input_LV(:, k) = xi_LV_k;

    if k == 1
        w_new = zeros(nx,1);
        alpha_before = alpha_ini;
        v_before = v_ini;
    else
        w_new = State_RM_UQMPC(:, k) - A*State_RM_UQMPC(:, k-1) - B*Control_EV_UQMPC(:, k-1);
        alpha_before = alpha_k;
        v_before = v_k;
    end

    tic
    [s_k_uqmpc, u_EV_k_uqmpc, alpha_k, v_k, hs_uqmpc, S_hat, W_hat] = UQRobustMPC.solve(State_RM_UQMPC(:, k), alpha_before, v_before, w_new);
    toc
    Time_UQMPC(k) = toc;

    tic
    [u_EV_k_scmpc] = SCMPC.solve(State_RM_SCMPC(:, k));
    toc
    Time_SCMPC(k) = toc;

    [xi_EV_k, ~] = Car.EVModeling(zeros(nx, 1), zeros(nu, 1)); % noise of EV

    x_EV_k_next_uqmpc = A*State_EV_UQMPC(:, k) + B*u_EV_k_uqmpc + xi_EV_k;
    x_EV_k_next_scmpc = A*State_EV_SCMPC(:, k) + B*u_EV_k_scmpc + xi_EV_k;

    Noise_EV(:, k)   = xi_EV_k;

    Control_EV_UQMPC(:, k)   = u_EV_k_uqmpc;
    State_EV_UQMPC(:, k + 1) = x_EV_k_next_uqmpc;
    State_RM_UQMPC(:, k + 1) = x_EV_k_next_uqmpc - x_LV_k_next - x_des;
    Nominal_RM_State_UQMPC(:, k) = s_k_uqmpc;
    Hs_UQMPC(:, k) = hs_uqmpc;
    S_Hat{k} = S_hat;
    W_Hat{k} = W_hat;

    Alpha(:, k + 1) = alpha_k;
    V{k + 1} = v_k;
    % Store samples as columns for any nx
    samples = [samples, w_new];

    Control_EV_SCMPC(:, k) = u_EV_k_scmpc;
    State_EV_SCMPC(:, k + 1) = x_EV_k_next_scmpc;
    State_RM_SCMPC(:, k + 1) = x_EV_k_next_scmpc - x_LV_k_next - x_des;
end

fprintf('Mean Excecution Time of UQ-RMPC is %.4f s.\n', sum(Time_UQMPC)/K_N);
fprintf('Max. Excecution Time of UQ-RMPC is %.4f s.\n', max(Time_UQMPC));

fprintf('Mean Excecution Time of SCMPC is %.4f. s\n', sum(Time_SCMPC)/K_N);
fprintf('Max. Excecution Time of SCMPC is %.4f s\n', max(Time_SCMPC));
%%
%{
figure(1)
if nx == 2
    plot(State_RM_UQMPC(1, :), State_RM_UQMPC(2, :), 'b', 'linewidth', 2)
    hold on
    plot(State_RM_SCMPC(1, :), State_RM_SCMPC(2, :), 'r--', 'linewidth', 1.5)
    xlabel('$x_{k,1}$ [${\rm m}$]', 'Interpreter','latex');
    ylabel('$x_{k,2}$ [${\rm m/s}$]', 'Interpreter','latex');
elseif nx == 3
    plot3(State_RM_UQMPC(1, :), State_RM_UQMPC(2, :), State_RM_UQMPC(3, :), 'b', 'linewidth', 2)
    hold on
plot3(State_RM_SCMPC(1, :), State_RM_SCMPC(2, :), State_RM_SCMPC(3, :), 'r--', 'linewidth', 1.5)
    xlabel('$x_{k,1}$', 'Interpreter','latex');
    ylabel('$x_{k,2}$', 'Interpreter','latex');
    zlabel('$x_{k,3}$', 'Interpreter','latex');
else
    % For nx>3, plot selected projections and warn user
    plot(State_RM_UQMPC(1, :), State_RM_UQMPC(2, :), 'b', 'linewidth', 2)
    hold on
    plot(State_RM_SCMPC(1, :), State_RM_SCMPC(2, :), 'r--', 'linewidth', 1.5)
    xlabel('$x_{k,1}$', 'Interpreter','latex');
    ylabel('$x_{k,2}$', 'Interpreter','latex');
    title('Warning: Only first two dimensions shown for nx > 3')
    % Optionally, plot more projections for higher dimensions
    if nx > 3
        figure;
        plot(State_RM_UQMPC(2, :), State_RM_UQMPC(3, :), 'b', 'linewidth', 2)
        hold on
        plot(State_RM_SCMPC(2, :), State_RM_SCMPC(3, :), 'r--', 'linewidth', 1.5)
        xlabel('$x_{k,2}$', 'Interpreter','latex');
        ylabel('$x_{k,3}$', 'Interpreter','latex');
        title('Projection: 2nd vs 3rd state component')
        % Additional projections for nx > 4
        if nx > 4
            figure;
            plot(State_RM_UQMPC(1, :), State_RM_UQMPC(3, :), 'b', 'linewidth', 2)
            hold on
            plot(State_RM_SCMPC(1, :), State_RM_SCMPC(3, :), 'r--', 'linewidth', 1.5)
            xlabel('$x_{k,1}$', 'Interpreter','latex');
            ylabel('$x_{k,3}$', 'Interpreter','latex');
            title('Projection: 1st vs 3rd state component')
        end
        % Optionally, show a warning for even higher nx
        if nx > 5
            warning('State dimension nx > 5. Only a few projections are shown. Consider dimensionality reduction for visualization.');
        end
    end
end
LE = legend('UQ-RMPC', 'SCMPC',  'Interpreter','latex', 'Location','best');
set(LE, 'Fontsize', 12);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 10 10]);
set(gcf, 'PaperSize', [16 7]);


figure(2)
if nu == 1
    plot(0:T:T*(K_N-1), Control_EV_UQMPC, 'b', 'linewidth', 2)
    hold on
    plot(0:T:T*(K_N-1), Control_EV_SCMPC, 'r', 'linewidth', 1.5)
    hold on
    plot(0:T:T*(K_N-1), 2*ones(1, K_N), 'k:', 'linewidth', 2)
    hold on
    plot(0:T:T*(K_N-1), -2*ones(1, K_N), 'k:', 'linewidth', 2)
    ylim([-3, 3]);
    ylabel('${\rm Long. \\ Acc.}$ [${\rm m/s}$]', 'Interpreter','latex');
else
    % For nu > 1, plot each input channel
    colors = lines(nu);
    for i = 1:nu
        plot(0:T:T*(K_N-1), Control_EV_UQMPC(i,:), 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', ['UQ-RMPC u_', num2str(i)]);
        hold on
        plot(0:T:T*(K_N-1), Control_EV_SCMPC(i,:), '--', 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', ['SCMPC u_', num2str(i)]);
    end
    % Optionally, plot bounds if known (here for [-2,2])
    plot(0:T:T*(K_N-1), 2*ones(1, K_N), 'k:', 'linewidth', 2, 'DisplayName', 'Upper Bound')
    plot(0:T:T*(K_N-1), -2*ones(1, K_N), 'k:', 'linewidth', 2, 'DisplayName', 'Lower Bound')
    ylim([-3, 3]);
    ylabel('${\rm Input}$', 'Interpreter','latex');
end
LE = legend('show', 'Interpreter','latex', 'Location','best');
xlabel('${\rm Time}$ [${\rm s}$]', 'Interpreter','latex');
set(LE, 'Fontsize', 12);
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 17 5]);
set(gcf, 'PaperSize', [16 7]);
%}