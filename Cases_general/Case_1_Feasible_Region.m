% This computes the feasible region of UQ-RMPC, corresponding to Fig. 2
clc
clear
close all
addpath('../Functions_general');
load('parameters_3x3.mat');
opts_ini_set = parameters.opts_ini_set;
W = opts_ini_set.W;
%%
N_sam_ini = [50; 500; 2000; 10000]; % change |I_0^w|
Alpha_ini = ones(length(N_sam_ini), 1);



Phi = parameters.opts_RMPC.Phi;
S = parameters.opts_RMPC.S;
Nu = parameters.opts_RMPC.Nu;

FeasibleRegion = ComputeFeasibleRegion(parameters.opts_feasible_region);

%%% CHANGE!
%V_ini = cell(length(N_sam_ini), 1); % Store v_ini as cell array (can be Polyhedron or vector)
%V_ini = ones(length(N_sam_ini), 1);
V_ini = cell(length(N_sam_ini), 1);


Samples_ini = cell(length(N_sam_ini), 1);
W_Hat_ini = cell(length(N_sam_ini), 1);
for k = 1:1:length(N_sam_ini)
    opts_ini_set.N_pre_sam = N_sam_ini(k);
    IniSet = InitialSetComputation(opts_ini_set); 
    [alpha_ini, v_ini, samples] = IniSet.solve(); % return \alpha_0, v_0, and samples
    Alpha_ini(k) = alpha_ini;
    V_ini{k} = v_ini;
    Samples_ini{k} = samples;
    W_Hat_ini{k} = (1 - alpha_ini)*v_ini + alpha_ini*W;
    % Compute S_hat_ini and F_N_Hat for each k
    S_hat_ini = alpha_ini*S + (1 - alpha_ini)*(inv(eye(size(Phi)) - Phi)*v_ini);
    hs_UQMPC = FeasibleRegion.ComFeasibleRegion_UQMPC_Hs_ini(S_hat_ini);
    [~, F_N_Hat] = FeasibleRegion.ComFeasibleRegion_UQMPC(S_hat_ini, hs_UQMPC);
    F_N_Hat_ini{k} = F_N_Hat;
end

Results_1.N_sam_ini = N_sam_ini;
Results_1.Alpha_ini = Alpha_ini;
Results_1.V_ini = V_ini;
Results_1.Samples_ini = Samples_ini;
Results_1.W_Hat_ini = W_Hat_ini;
Results_1.W = W;
Results_1.W_hat_opt = parameters.W_hat_opt;
Results_1.W_true = parameters.W_true;
Results_1.F_N_True = parameters.F_N_True;
Results_1.F_N_Hat_Opt = parameters.F_N_Hat_Opt;
Results_1.F_N_RMPC = parameters.F_N_RMPC;
Results_1.F_N_Hat_ini = F_N_Hat_ini;
save('Results_1_3x3.mat', 'Results_1')
