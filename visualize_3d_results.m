% Visualize 3D Results from UQ-RMPC
% This script loads and visualizes 3D data for uncertainty sets and feasible regions.

clc
clear
close all

% --- Load Data ---
% Make sure the .mat files are in the MATLAB path or current directory
disp('Loading 3D results data...');
try
    % Load the specific 3D results file
    load('Results_1_3x3.mat');
    % The corresponding parameters file should have a matching name
    load('parameters_3x3.mat');
catch ME
    error('Could not load the required .mat files. Make sure "Results_1_3x3.mat" and "parameters_3x3.mat" are in the current path. Error: %s', ME.message);
end
disp('Data loaded successfully.');

% Extract variables from loaded data
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
F_N_Hat_ini = Results_1.F_N_Hat_ini;

Phi = parameters.opts_RMPC.Phi;
S = parameters.opts_RMPC.S;
Nu = parameters.opts_RMPC.Nu;

% --- Plotting Setup ---
FeasibleRegion = ComputeFeasibleRegion(parameters.opts_feasible_region);
c_w    = [255, 63, 164]/255;
c_hat  = [255, 163, 60]/255;
c_true = [177, 94, 255]/255;
c_opt  = [61, 48, 162]/255;
c5     = [88, 163, 153]/255;
c6 = [247, 237, 42]/255; % yellow
%{%
% --- Plot 1: Uncertainty Sets (3D) ---
disp('Generating 3D plots for uncertainty sets...');
for i = 1:length(N_sam_ini)
    figure('Name', sprintf('Uncertainty Sets 3D (N_sam=$%d)', N_sam_ini(i)));
    hold on;    

    % Plot samples in 3D
    samples = Samples_ini{i};
    if size(samples, 1) == 3
        h5 = plot3(samples(1, :), samples(2, :), samples(3, :), 'color', c5, 'marker', '.','markersize', 10, 'LineStyle','none');
    else
        warning('Sample data is not 3D, skipping sample plot.');
    end

    % Plot polyhedra (MPT's plot command handles 3D automatically)
    h1 = plot(W, 'wire', 1, 'edgecolor', c_w, 'linewidth', 3.0, 'linestyle', '-');
    h2 = plot(W_true, 'wire', 1, 'edgecolor', c_true, 'linewidth', 3.0, 'linestyle', '-');
    h4 = plot(W_Hat_ini{i}, 'wire', 1, 'edgecolor', c_hat, 'linewidth', 3.0);

    % Plot bounding cuboid for W_hat_opt instead of full polytope
    if isa(W_hat_opt, 'Polyhedron') && isprop(W_hat_opt, 'V') && size(W_hat_opt.V,2)==3 && size(W_hat_opt.V,1)>=3
        V = W_hat_opt.V;
        minV = min(V,[],1);
        maxV = max(V,[],1);
        % Create cuboid vertices
        cuboidVerts = [minV;
                       minV(1), minV(2), maxV(3);
                       minV(1), maxV(2), minV(3);
                       minV(1), maxV(2), maxV(3);
                       maxV(1), minV(2), minV(3);
                       maxV(1), minV(2), maxV(3);
                       maxV(1), maxV(2), minV(3);
                       maxV(1), maxV(2), maxV(3)];
        % Define cuboid faces
        faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
        % Plot cuboid as semi-transparent patch (wireframe)
        h3 = patch('Vertices', cuboidVerts, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', c_opt, 'LineWidth', 3, 'LineStyle', '-', 'DisplayName', '$\widehat{W}_{\rm opt}$ (cuboid)');
    else
        h3 = plot(W_hat_opt, 'wire', 1, 'edgecolor', c_opt, 'linewidth', 3.0, 'linestyle', '-');
    end

    

    hold off;
    box on;
    grid on;
    view(3); % Set 3D view
    axis equal;
    
    xlabel('$w_{k, 1}$', 'Interpreter','latex');
    ylabel('$w_{k, 2}$', 'Interpreter','latex');
    zlabel('$w_{k, 3}$', 'Interpreter','latex');
    title(['Uncertainty Sets ($N_{sam}=' num2str(N_sam_ini(i)) '$)'], 'Interpreter', 'latex');    
    set(gca,'Linewidth',1.0,'GridAlpha',0.5);
    set(gca,'FontName','Times New Roman','FontSize',15);
    legend([h1 h2 h3 h4 h5], {'$W$', '$W_{\rm true}$', '$\widehat{W}_{\rm opt}$', '$\widehat{W}_0^*$','Disturbance samples$'}, 'Interpreter', 'latex', 'Location', 'best');
    
    savename = sprintf('Fig_3D_W_hat_vs_W_Samples_%d.pdf', N_sam_ini(i));
    %exportgraphics(gcf, savename, 'ContentType', 'vector');
end


% --- Additional Figure: Overlay W_Hat_ini{1} and W_Hat_ini{2} on top of W, W_true, W_hat_opt ---
figure('Name', 'Uncertainty Sets 3D with Initial Sets Overlay');
hold on;
samples = Samples_ini{4};
h6 = plot3(samples(1, :), samples(2, :), samples(3, :), 'color', c_hat, 'marker', '.','markersize', 10, 'LineStyle','none');
samples = Samples_ini{1};
h7 = plot3(samples(1, :), samples(2, :), samples(3, :), 'color', c5, 'marker', 'x','markersize', 15, 'LineStyle','none', 'linewidth', 2);
% Plot main sets as in Plot 1 (using cuboid for W_hat_opt)
h1 = plot(W, 'wire', 1, 'edgecolor', c_w, 'linewidth', 3.0, 'linestyle', '-');
h2 = plot(W_true, 'wire', 1, 'edgecolor', c_true, 'linewidth', 3.0, 'linestyle', '-');

if isa(W_hat_opt, 'Polyhedron') && isprop(W_hat_opt, 'V') && size(W_hat_opt.V,2)==3 && size(W_hat_opt.V,1)>=3
    V = W_hat_opt.V;
    minV = min(V,[],1);
    maxV = max(V,[],1);
    cuboidVerts = [minV;
                   minV(1), minV(2), maxV(3);
                   minV(1), maxV(2), minV(3);
                   minV(1), maxV(2), maxV(3);
                   maxV(1), minV(2), minV(3);
                   maxV(1), minV(2), maxV(3);
                   maxV(1), maxV(2), minV(3);
                   maxV(1), maxV(2), maxV(3)];
    faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
    h3 = patch('Vertices', cuboidVerts, 'Faces', faces, 'FaceColor', 'none', 'EdgeColor', c_opt, 'LineWidth', 3, 'LineStyle', '-', 'DisplayName', '$\\widehat{W}_{\\rm opt}$ (cuboid)');
else
    h3 = plot(W_hat_opt, 'wire', 1, 'edgecolor', c_opt, 'linewidth', 3.0, 'linestyle', '-');
end
% Overlay W_Hat_ini{1} and W_Hat_ini{2} with different styles
if length(W_Hat_ini) >= 1
    h4 = plot(W_Hat_ini{1}, 'wire', 1, 'edgecolor', c5, 'linewidth', 3.0, 'linestyle', '-');
else
    h4 = [];
end
if length(W_Hat_ini) >= 2
    h5 = plot(W_Hat_ini{4}, 'wire', 1, 'edgecolor', c_hat, 'linewidth', 3.0, 'linestyle', '-');
else
    h5 = [];
end

hold off;
box on;
grid on;
view(3);
axis equal;
xlabel('$w_{k, 1}$', 'Interpreter','latex');
ylabel('$w_{k, 2}$', 'Interpreter','latex');
zlabel('$w_{k, 3}$', 'Interpreter','latex');
title('Uncertainty Sets with Initial Sets Overlay', 'Interpreter', 'latex');
set(gca,'Linewidth',1.0,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',15);
legend([h1 h2 h3 h4 h5 h7 h6], {'$W$', '$W_{\rm true}$', '$\widehat{W}_{\rm opt}$', '$\widehat{W}_{0,50}^*$', '$\widehat{W}_{0,2000}^*$', '$|\mathcal{I}_0^w|=50$','$|\mathcal{I}_0^w|=2000$'}, 'Interpreter', 'latex', 'Location', 'best');
%}


% --- Plot 2: Feasible Regions (3D) ---
disp('Generating 3D plots for feasible regions...');
%{
% --- ORIGINAL CODE: Re-calculates feasible regions, which is very slow ---
for i = 1:length(N_sam_ini)
    figure('Name', sprintf('Feasible Regions 3D (N_sam=%d)', N_sam_ini(i)));
    hold on;

    % Calculate initial feasible set for this iteration
    S_hat_ini = Alpha_ini(i)*S + (1 - Alpha_ini(i))*(inv(eye(size(Phi)) - Phi)*V_ini(:, i));
    hs_UQMPC = FeasibleRegion.ComFeasibleRegion_UQMPC_Hs_ini(S_hat_ini);
    [~, F_N_Hat] = FeasibleRegion.ComFeasibleRegion_UQMPC(S_hat_ini, hs_UQMPC);
    
    % Plot polyhedra
    plot(F_N_True, 'wire', 1, 'edgecolor', c_true, 'linewidth', 1.5, 'linestyle', '--');
    plot(F_N_Hat, 'wire', 1, 'edgecolor', c_hat, 'linewidth', 1.5);
    plot(F_N_RMPC, 'wire', 1, 'edgecolor', c_w, 'linewidth', 1.5, 'linestyle', '-');
    plot(F_N_Hat_Opt, 'wire', 1, 'edgecolor', c_opt, 'linewidth', 2.0, 'linestyle', '-.');
    
    hold off;
    box on;
    grid on;
    view(3); % Set 3D view
    axis tight;

    xlabel('$x_{0,1}$', 'Interpreter','latex');
    ylabel('$x_{0,2}$', 'Interpreter','latex');
    zlabel('$x_{0,3}$', 'Interpreter','latex');
    title(sprintf('Feasible Regions (N_{sam}=%d)', N_sam_ini(i)), 'Interpreter', 'none');

    set(gca,'Linewidth',1.5,'GridAlpha',0.5);
    set(gca,'FontName','Times New Roman','FontSize',15);
    
    savename = sprintf('Fig_3D_F_N_hat_vs_F_N_%d.pdf', N_sam_ini(i));
    %exportgraphics(gcf, savename, 'ContentType', 'vector');
end
%}

% --- NEW CODE: Plots only the final, pre-computed feasible regions ---
figure('Name', 'Feasible Regions 3D (Final)');
hold on;

% --- DEBUG STEP: Check content of feasible region polyhedra ---
disp('--- Debugging Feasible Regions ---');
fprintf('Variable: F_N_True\n');
disp(F_N_True);
if isa(F_N_True, 'Polyhedron') && F_N_True.Dim > 0
    fprintf('Is empty: %d, Is bounded: %d, Vertices: %d, Constraints: %d\n\n', ...
        F_N_True.isEmptySet, F_N_True.isBounded, size(F_N_True.V, 1), size(F_N_True.A, 1));
end

fprintf('Variable: F_N_RMPC\n');
disp(F_N_RMPC);
if isa(F_N_RMPC, 'Polyhedron') && F_N_RMPC.Dim > 0
    fprintf('Is empty: %d, Is bounded: %d, Vertices: %d, Constraints: %d\n\n', ...
        F_N_RMPC.isEmptySet, F_N_RMPC.isBounded, size(F_N_RMPC.V, 1), size(F_N_RMPC.A, 1));
end

fprintf('Variable: F_N_Hat_Opt\n');
disp(F_N_Hat_Opt);
if isa(F_N_Hat_Opt, 'Polyhedron') && F_N_Hat_Opt.Dim > 0
    fprintf('Is empty: %d, Is bounded: %d, Vertices: %d, Constraints: %d\n\n', ...
        F_N_Hat_Opt.isEmptySet, F_N_Hat_Opt.isBounded, size(F_N_Hat_Opt.V, 1), size(F_N_Hat_Opt.A, 1));
end
disp('--- End Debugging ---');

% --- Polyhedron vertex reduction utility (toolbox-free) ---
% This function reduces the number of vertices by random sampling.
% Usage: Vred = reduce_vertices(V, max_vertices)
% V: Nx3 matrix of vertices
% max_vertices: desired maximum number of vertices
function Vred = reduce_vertices(V, max_vertices)
    if size(V,1) <= max_vertices
        Vred = V;
        return;
    end
    % Randomly sample max_vertices from V
    idx = randperm(size(V,1), max_vertices);
    Vred = V(idx, :);
end

% --- Efficient Polyhedron Visualization with Bounding Wireframe ---
% (Original code is commented out below for reference)
obbox_colors = {c_true, c_w, c_opt, c_hat, c5};
obbox_names = {'$F_N^{\rm true}$', '$F_N^{\rm RMPC}$', '$\widehat{F}_N^{\rm opt}$', '$\widehat{F}^{\rm 1}$', '$\widehat{F}^{\rm 4}$'};

max_vis_vertices = 35; % You can adjust this value for speed/quality tradeoff
h1 = [];
h2 = [];
h3 = [];
h4 = [];
h5 = [];
for k = 1:3
    switch k
        case 1, P = F_N_True; col = c_true; lw = 2; ls = '-.'; name = 'F_N_True';
        case 2, P = F_N_RMPC; col = c_w; lw = 2; ls = '-'; name = 'F_N_RMPC';
        case 3, P = F_N_Hat_ini{4}; col = obbox_colors{4}; lw = 2; name = obbox_names{4};
        case 4, P = F_N_Hat_Opt; col = c_opt; lw = 2; ls = '-'; name = 'F_N_Hat_Opt';
        case 5, P = F_N_Hat_ini{1}; col = obbox_colors{5}; lw = 1; name = obbox_names{5};
    end
    if isa(P, 'Polyhedron') && P.Dim == 3 && isprop(P, 'V') && size(P.V,1) > max_vis_vertices
        disp(['Reducing vertices for ' name ' from ' num2str(size(P.V,1)) ' to ' num2str(max_vis_vertices)]);
        Vred = reduce_vertices(P.V, max_vis_vertices);
        % Create a new Polyhedron with reduced vertices (convex hull)
        Pred = Polyhedron('V', Vred);
        % Plot only the wireframe of the bounding polyhedron
        h = plot(Pred, 'wire', 1, 'edgecolor', col, 'linewidth', lw, 'linestyle', ls);
    elseif isa(P, 'Polyhedron') && P.Dim == 3
        disp(['Plotting ' name ' as bounding wireframe (no reduction needed).']);
        h = plot(P, 'wire', 1, 'edgecolor', col, 'linewidth', lw, 'linestyle', ls);
    else
        disp(['Plotting ' name ' as wireframe or skipping (not 3D Polyhedron).']);
        h = plot(P, 'wire', 1, 'edgecolor', col, 'linewidth', lw, 'linestyle', ls);
    end
    switch k
        case 1, h1 = h;
        case 2, h2 = h;
        case 3, h3 = h;
        case 4, h4 = h;
        case 5, h5 = h;
    end
end

legend([h1 h2 h3 h4 h5], {'$F_N^{\rm true}$', '$F_N^{\rm RMPC}$', '$\widehat{F}_{0,2000}^{\rm *}$', '$\widehat{F}_N^{\rm opt}$', '$\widehat{F}_{0,50}^{\rm *}$'}, 'Interpreter', 'latex', 'Location', 'best');

%{
% --- Save 2D cross-sections of 3D polyhedra at x1=0, x2=0, x3=0 planes ---
disp('Saving 2D cross-sections of feasible regions at x1=0, x2=0, x3=0...');
section_planes = [1 2 3]; % x1, x2, x3
section_names = {'x1', 'x2', 'x3'};
poly_names = {'F_N_True', 'F_N_RMPC', 'F_N_Hat_Opt'};
poly_list = {F_N_True, F_N_RMPC, F_N_Hat_Opt};
poly_colors = {c_true, c_w, c_opt};

for pidx = 1:length(poly_list)
    P = poly_list{pidx};
    pname = poly_names{pidx};
    pcolor = poly_colors{pidx};
    if isa(P, 'Polyhedron') && isprop(P, 'V') && size(P.V,2)==3 && size(P.V,1)>=3
        V = P.V;
        edges = nchoosek(1:size(V,1), 2);
        for sidx = 1:length(section_planes)
            dim = section_planes(sidx);
            val = 0; % plane value
            section_pts = [];
            for eidx = 1:size(edges,1)
                v1 = V(edges(eidx,1),:);
                v2 = V(edges(eidx,2),:);
                if (v1(dim) - val) * (v2(dim) - val) < 0
                    t = (val - v1(dim)) / (v2(dim) - v1(dim));
                    pt = v1 + t * (v2 - v1);
                    section_pts = [section_pts; pt];
                end
            end
            % Plot section in other two dimensions
            figure('Name', sprintf('Section of %s at %s=0', pname, section_names{sidx}));
            hold on;
            if dim == 1
                plot(section_pts(:,2), section_pts(:,3), '.', 'Color', pcolor, 'MarkerSize', 15);
                xlabel('x_2'); ylabel('x_3');
            elseif dim == 2
                plot(section_pts(:,1), section_pts(:,3), '.', 'Color', pcolor, 'MarkerSize', 15);
                xlabel('x_1'); ylabel('x_3');
            else
                plot(section_pts(:,1), section_pts(:,2), '.', 'Color', pcolor, 'MarkerSize', 15);
                xlabel('x_1'); ylabel('x_2');
            end
            title(sprintf('Section of %s at %s = 0', pname, section_names{sidx}));
            axis equal; grid on;
            set(gca,'FontName','Times New Roman','FontSize',15);
            hold off;
            fname = sprintf('Section_%s_%s0.png', pname, section_names{sidx});
            exportgraphics(gcf, fname, 'Resolution', 300);
            disp(['Saved section: ' fname]);
        end
    else
        disp(['Skipping section for ' pname ': not a 3D Polyhedron']);
    end
end
%}

% (Legend unchanged)

% --- Original ellipsoid/polyhedron plotting code is commented out below for reference ---
% --- Ellipsoid-only Visualization ---
% disp('Ellipsoid-only visualization: fitting and plotting ellipsoids for each region.');
% h1 = [];
% h2 = [];
% h3 = [];
% for k = 1:3
%     switch k
%         case 1, P = F_N_True; col = c_true; alp = 0.10; name = 'F_N_True';
%         case 2, P = F_N_RMPC; col = c_w; alp = 0.13; name = 'F_N_RMPC';
%         case 3, P = F_N_Hat_Opt; col = c_opt; alp = 0.16; name = 'F_N_Hat_Opt';
%     end
%     if isa(P, 'Polyhedron') && isprop(P, 'V') && size(P.V,2)==3 && size(P.V,1)>=3
%         disp(['Fitting ellipsoid for ' name '...']);
%         try
%             [A, c] = minVolEllipse(P.V', 1e-3); % A: shape matrix, c: center
%             h = plot_ellipsoid(A, c, col, alp);
%         catch ME
%             warning(['Ellipsoid fit failed for ' name ': ' ME.message]);
%             h = [];
%         end
%     else
%         warning(['Skipping ellipsoid for ' name ': not enough vertices or not 3D.']);
%         h = [];
%     end
%     switch k
%         case 1, h1 = h;
%         case 2, h2 = h;
%         case 3, h3 = h;
%     end
% end
% legend([h1 h2 h3], {'$F_N^{\\rm true}$', '$F_N^{\\rm RMPC}$', '$\\widehat{F}_N^{\\rm opt}$'}, 'Interpreter', 'latex', 'Location', 'best');
% (Legend unchanged)

% --- Helper function for plotting ellipsoids ---
function h = plot_ellipsoid(A, c, color, alpha)
    [X,Y,Z] = ellipsoid(0,0,0,1,1,1,30);
    XYZ = [X(:) Y(:) Z(:)]';
    [U,S,~] = svd(A);
    radii = 1./sqrt(diag(S));
    ell = U*diag(radii)*XYZ + c;
    X2 = reshape(ell(1,:), size(X));
    Y2 = reshape(ell(2,:), size(Y));
    Z2 = reshape(ell(3,:), size(Z));
    h = surf(X2, Y2, Z2, 'FaceColor', color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end

hold off;
box on;
grid on;
view(3); % Set 3D view
axis tight;

xlabel('$x_{0,1}$', 'Interpreter','latex');
ylabel('$x_{0,2}$', 'Interpreter','latex');
zlabel('$x_{0,3}$', 'Interpreter','latex');
title('Final Feasible Regions', 'Interpreter', 'none');

set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',15);
% legend([h1 h2 h3], {'$F_N^{\rm true}$', '$F_N^{\rm RMPC}$', '$\widehat{F}_N^{\rm opt}$'}, 'Interpreter', 'latex', 'Location', 'best');
% (Legend unchanged)

savename = 'Fig_3D_Feasible_Regions_Final.pdf';
%exportgraphics(gcf, savename, 'ContentType', 'vector');


% --- Bounding Cuboid Visualization ---
figure('Name', 'Bounding Cuboids for Feasible Regions');
hold on;
cuboid_colors = {c_true, c_w, c_opt};
cuboid_names = {'$F_N^{\rm true}$', '$F_N^{\rm RMPC}$', '$\widehat{F}_N^{\rm opt}$'};
h_cub = [];
for k = 1:3
    switch k
        case 1, P = F_N_True; col = cuboid_colors{1}; name = cuboid_names{1};
        case 2, P = F_N_RMPC; col = cuboid_colors{2}; name = cuboid_names{2};
        case 3, P = F_N_Hat_Opt; col = cuboid_colors{3}; name = cuboid_names{3};
    end
    if isa(P, 'Polyhedron') && isprop(P, 'V') && size(P.V,2)==3 && size(P.V,1)>=3
        V = P.V;
        minV = min(V,[],1);
        maxV = max(V,[],1);
        % Create cuboid vertices
        cuboidVerts = [minV;
                       minV(1), minV(2), maxV(3);
                       minV(1), maxV(2), minV(3);
                       minV(1), maxV(2), maxV(3);
                       maxV(1), minV(2), minV(3);
                       maxV(1), minV(2), maxV(3);
                       maxV(1), maxV(2), minV(3);
                       maxV(1), maxV(2), maxV(3)];
        % Define cuboid faces
        faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
        % Plot cuboid as semi-transparent patch
        h = patch('Vertices', cuboidVerts, 'Faces', faces, 'FaceColor', col, 'FaceAlpha', 0.00, 'EdgeColor', col, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', name);
        h_cub = [h_cub h];
    end
end
hold off;
grid on;
box on;
view(3);
axis tight;
xlabel('$x_{0,1}$', 'Interpreter','latex');
ylabel('$x_{0,2}$', 'Interpreter','latex');
zlabel('$x_{0,3}$', 'Interpreter','latex');
title('Bounding Cuboids for Final Feasible Regions', 'Interpreter', 'none');
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',15);
legend(h_cub, cuboid_names, 'Interpreter', 'latex', 'Location', 'best');


% --- Oblique (minimum volume) bounding box visualization ---
figure('Name', 'Oblique Bounding Boxes for Feasible Regions');
hold on;
obbox_colors = {c_true, c_w, c_opt, c_hat, c5};
obbox_names = {'$F_N^{\rm true}$', '$F_N^{\rm RMPC}$', '$\widehat{F}_N^{\rm opt}$', '$\widehat{F}_{0,2000}^{\rm*}$', '$\widehat{F}__{0,50}^{\rm *}$'};
h_obbox = [];
for k = 1:5
    switch k
        case 1, P = F_N_True; col = obbox_colors{1}; name = obbox_names{1};
        case 2, P = F_N_RMPC; col = obbox_colors{2}; name = obbox_names{2};
        case 3, P = F_N_Hat_Opt; col = obbox_colors{3}; name = obbox_names{3};
        case 4, P = F_N_Hat_ini{1}; col = obbox_colors{4}; name = obbox_names{4};
        case 5, P = F_N_Hat_ini{4}; col = obbox_colors{5}; name = obbox_names{5};

    end
    if isa(P, 'Polyhedron') && isprop(P, 'V') && size(P.V,2)==3 && size(P.V,1)>=3
        V = P.V;
        % PCA for oblique bounding box
        [coeff,score,~] = pca(V);
        Vp = V*coeff; % rotate to principal axes
        %Vp=V;
        minVp = min(Vp,[],1);
        maxVp = max(Vp,[],1);
        % Vertices of the box in PCA space
        boxVerts = [minVp;
                    minVp(1), minVp(2), maxVp(3);
                    minVp(1), maxVp(2), minVp(3);
                    minVp(1), maxVp(2), maxVp(3);
                    maxVp(1), minVp(2), minVp(3);
                    maxVp(1), minVp(2), maxVp(3);
                    maxVp(1), maxVp(2), minVp(3);
                    maxVp(1), maxVp(2), maxVp(3)];
        % Transform back to original space
        boxVertsWorld = boxVerts*coeff';
        %boxVertsWorld = boxVerts;
        faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
        h = patch('Vertices', boxVertsWorld, 'Faces', faces, 'FaceColor', col, 'FaceAlpha', 0.00, 'EdgeColor', col, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', name);
        h_obbox = [h_obbox h];
    end
end
hold off;
grid on;
box on;
view(3);
axis tight;
xlabel('$x_{0,1}$', 'Interpreter','latex');
ylabel('$x_{0,2}$', 'Interpreter','latex');
zlabel('$x_{0,3}$', 'Interpreter','latex');
title('Oblique Bounding Boxes for Final Feasible Regions', 'Interpreter', 'none');
set(gca,'Linewidth',1.5,'GridAlpha',0.5);
set(gca,'FontName','Times New Roman','FontSize',15);
legend(h_obbox, obbox_names, 'Interpreter', 'latex', 'Location', 'best');

% --- Print volumes of bounding boxes for feasible regions ---
for k = 1:3
    switch k
        case 1, P = F_N_True; name = 'F_N_True';
        case 2, P = F_N_RMPC; name = 'F_N_RMPC';
        case 3, P = F_N_Hat_Opt; name = 'F_N_Hat_Opt';
    end
    if isa(P, 'Polyhedron') && isprop(P, 'V') && size(P.V,2)==3 && size(P.V,1)>=3
        V = P.V;
        minV = min(V,[],1);
        maxV = max(V,[],1);
        box_vol = prod(maxV - minV);
        fprintf('Bounding box volume for %s: %.4f\n', name, box_vol);
    else
        fprintf('Bounding box volume for %s: N/A (not 3D Polyhedron)\n', name);
    end
end

disp('All visualizations generated and saved as PDF files.');

% --- minVolEllipse utility (Khachiyan Algorithm, public domain) ---
% Source: https://www.mathworks.com/matlabcentral/fileexchange/9542-minimum-volume-enclosing-ellipsoid
% Usage: [A, c] = minVolEllipse(P, tolerance)
% P: 3xN matrix of points (each column is a point)
% tolerance: e.g. 1e-3
% Returns: A (shape matrix), c (center)
function [A, c] = minVolEllipse(P, tolerance)
    [d, N] = size(P);
    Q = [P; ones(1,N)];
    count = 1;
    err = 1;
    u = (1/N) * ones(N,1); % initial weights
    while err > tolerance
        X = Q * diag(u) * Q';
        M = diag(Q' / X * Q);
        [~, j] = max(M);
        step_size = (M(j) - d - 1)/((d+1)*(M(j)-1));
        new_u = (1-step_size)*u; new_u(j) = new_u(j) + step_size;
        err = norm(new_u-u);
        u = new_u;
        count = count + 1;
        if count > 1e4
            warning('minVolEllipse: Max iterations reached');
            break;
        end
    end
    c = P * u;
    A = inv(P*diag(u)*P' - c*c')/d;
end

% --- End minVolEllipse utility ---


%}