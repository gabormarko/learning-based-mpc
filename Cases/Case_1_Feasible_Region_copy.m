%% This computes the feasible region of UQ-RMPC, corresponding to Fig. 2
clc
clear
close all
load('parameters_NEW.mat');
opts_ini_set = parameters.opts_ini_set;
W = opts_ini_set.W;
%%
N_sam_ini = [20; 200; 20000]; % change |I_0^w|
Alpha_ini = ones(length(N_sam_ini), 1);

%%% CHANGE!
V_ini = ones(parameters.opts_ini_set.nx, length(N_sam_ini));

Samples_ini = cell(length(N_sam_ini), 1);
W_Hat_ini = cell(length(N_sam_ini), 1);
for k = 1:1:length(N_sam_ini)
    opts_ini_set.N_pre_sam = N_sam_ini(k);
    IniSet = InitialSetComputation(opts_ini_set); 
    [alpha_ini, v_ini, samples] = IniSet.solve(); % return \alpha_0, v_0, and samples
    Alpha_ini(k) = alpha_ini;
    V_ini(:, k) = v_ini;
    fprintf('N_sam_ini = %d, alpha_ini = %.4f, v_ini = %.4f\n', N_sam_ini(k), alpha_ini, v_ini);
    Samples_ini{k} = samples;
    W_Hat_ini{k} = (1 - alpha_ini)*v_ini + alpha_ini*W;
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
save('Results_1_NEW.mat', 'Results_1')



% Load parameters
%clc
%clear
close all
load('Results_1_NEW.mat');
load('parameters_NEW.mat');
N_sam_ini = Results_1.N_sam_ini;
Samples_ini = Results_1.Samples_ini;
Alpha_ini = Results_1.Alpha_ini;
V_ini = Results_1.V_ini;
W_Hat_ini = Results_1.W_Hat_ini;
W = Results_1.W;
W_hat_opt = Results_1.W_hat_opt;
W_true = Results_1.W_true;
F_N_True = Results_1.F_N_True;
F_N_Hat_Opt = Results_1.F_N_Hat_Opt;
F_N_RMPC = Results_1.F_N_RMPC;

Phi = parameters.opts_RMPC.Phi;
S = parameters.opts_RMPC.S;
Nu = parameters.opts_RMPC.Nu;

FeasibleRegion = ComputeFeasibleRegion(parameters.opts_feasible_region);
c_w    = [219, 134, 31]/255;    %orange
c_hat2  = [87, 172, 242]/255;  %blue
c_hat1  = [181, 255, 254]/255;  % light blue
c_true = [31, 219, 65]/255;     % green
c_opt  = [29, 2, 163]/255;      % dark blue
c51     = [166, 31, 219]/255;
c52     = [255, 255, 89]/255;     
%%
close all
figure(1)
h1 = plot(W, 'wire', 1, 'edgecolor', c_w, 'linewidth', 2.5);
hold on
h6 = plot(W_Hat_ini{1}, 'wire', 1, 'edgecolor', c_hat1, 'linewidth', 2.5);
hold on
h7 = plot(W_Hat_ini{3}, 'wire', 1, 'edgecolor', c_hat2, 'linewidth', 2.5);
hold on
h2 = plot(W_true, 'wire', 1, 'edgecolor', c_true, 'linewidth', 2.5);
hold on
h3 = plot(W_hat_opt, 'wire', 1, 'edgecolor', c_opt, 'linewidth', 2.5);
hold on
h4 = plot(Samples_ini{1}(1, :),Samples_ini{1}(2, :), 'color', c51, 'marker', '.','markersize', 10, 'LineStyle','none');
box on
h5 = plot(Samples_ini{3}(1, :),Samples_ini{3}(2, :), 'color', c52, 'marker', '.','markersize', 10, 'LineStyle','none');
box on
grid off
xlim([-0.62, 0.62]);
ylim([-0.25, 0.25]);
xlabel('$w_{k, 1}$ [${\rm m}$]', 'Interpreter','latex');
ylabel('$w_{k, 2}$ [${\rm m/s}$]', 'Interpreter','latex');
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',15);
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'unit','centimeters','position',[5 5 10 10]);
set(gcf, 'PaperSize', [16 7]);
%legend([h1 h2 h3 h4 h5 h6 h7], {'$W$', '$W_{\rm true}$', '$\widehat{W}_{\rm opt}$', '\abs{\mathcal{I}_0^w}=20', '\abs{\mathcal{I}_0^w}=20000', '$\widehat{W}_{0,20}^*$', '$\widehat{W}_{0,20000}^*$'}, 'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 14);
savename = sprintf('NEW_Fig_W_hat_vs_W_Samples_%d.pdf', N_sam_ini(3));
exportgraphics(gcf, savename,'ContentType','vector');

%%
close all
for i = 1:length(N_sam_ini)
    figure(i)
    S_hat_ini = Alpha_ini(i)*S + (1 - Alpha_ini(i))*(inv(1 - Phi)*V_ini(:, i));
    hs_UQMPC = FeasibleRegion.ComFeasibleRegion_UQMPC_Hs_ini(S_hat_ini);
    [~, F_N_Hat] = FeasibleRegion.ComFeasibleRegion_UQMPC(S_hat_ini, hs_UQMPC);
    plot(F_N_True, 'wire', 1, 'edgecolor', c_true, 'linewidth', 2.5);
    hold on
    plot(F_N_Hat, 'wire', 1, 'edgecolor', c_hat, 'linewidth', 2.5);
    hold on
    plot(F_N_RMPC, 'wire', 1, 'edgecolor', c_w, 'linewidth', 2.5);
    hold on
    plot(F_N_Hat_Opt, 'wire', 1, 'edgecolor', c_opt, 'linewidth', 2.5);
    box on
    grid off
    xlim([-16, 16]);
    ylim([-8, 8]);
    xlabel('$x_{0,1}$ [${\rm m}$]', 'Interpreter','latex');
    ylabel('$x_{0, 2}$ [${\rm m/s}$]', 'Interpreter','latex');
    set(gca,'Linewidth',1.5,'GridAlpha',0.5);
    set(gca,'FontName','Times New Roman','FontSize',15);
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'unit','centimeters','position',[5 5 10 10]);
    set(gcf, 'PaperSize', [16 7]);
    savename = sprintf('NEW_Fig_F_N_hat_vs_F_N_%d.pdf', N_sam_ini(i));
    exportgraphics(gcf, savename,'ContentType','vector');
end