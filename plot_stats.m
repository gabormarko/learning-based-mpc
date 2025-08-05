% plot_stats.m
% Plot statistics for different system dimensions (nx)
% Place this script in the same folder as your stats_nxX.mat files

dims = [2, 3]%, 4]; % Add more if you run higher dimensions
all_stats = cell(length(dims),1);

for i = 1:length(dims)
    fname = ['stats_nx', num2str(dims(i)), '.mat'];
    if exist(fname, 'file')
        load(fname, 'stats');
        all_stats{i} = stats;
    else
        warning(['File not found: ', fname]);
        all_stats{i} = [];
    end
end

% Remove empty entries
valid = ~cellfun(@isempty, all_stats);
dims = dims(valid);
all_stats = all_stats(valid);

%% Plot computation times
figure;
hold on;
plot(dims, cellfun(@(s) s.t_mrpi_shat, all_stats), '-o', 'DisplayName', 'MRPI S\_hat\_opt');
plot(dims, cellfun(@(s) s.t_mrpi_strue, all_stats), '-o', 'DisplayName', 'MRPI S\_true');
plot(dims, cellfun(@(s) s.t_mrpi_s, all_stats), '-o', 'DisplayName', 'MRPI S (scaled)');
plot(dims, cellfun(@(s) s.t_feareg, all_stats), '-o', 'DisplayName', 'Feasible Region');
plot(dims, cellfun(@(s) s.t_sdp, all_stats), '-o', 'DisplayName', 'SDP for Ps');
plot(dims, cellfun(@(s) s.t_initialset, all_stats), '-o', 'DisplayName', 'InitialSetComputation');
plot(dims, cellfun(@(s) s.elapsedTime, all_stats), '-o', 'DisplayName', 'Total');
xlabel('System Dimension (nx)');
ylabel('Time (s)');
legend('Location','northwest');
title('Computation Time vs. System Dimension');
grid on;
hold off;

%% Plot number of vertices in S_hat_opt
figure;
plot(dims, cellfun(@(s) s.S_hat_opt_num_vertices, all_stats), '-o');
xlabel('System Dimension (nx)');
ylabel('Number of Vertices');
title('Vertices in S\_hat\_opt vs. Dimension');
grid on;

%% Plot number of constraints in S_hat_opt
figure;
plot(dims, cellfun(@(s) s.S_hat_opt_num_constraints, all_stats), '-o');
xlabel('System Dimension (nx)');
ylabel('Number of Constraints');
title('Constraints in S\_hat\_opt vs. Dimension');
grid on;

%% Plot feasible region sizes
figure;
plot(dims, cellfun(@(s) prod(s.F_N_RMPC_size), all_stats), '-o', 'DisplayName', 'F\_N\_RMPC size');
hold on;
plot(dims, cellfun(@(s) prod(s.F_N_True_size), all_stats), '-o', 'DisplayName', 'F\_N\_True size');
plot(dims, cellfun(@(s) prod(s.F_N_Hat_Opt_size), all_stats), '-o', 'DisplayName', 'F\_N\_Hat\_Opt size');
xlabel('System Dimension (nx)');
ylabel('Feasible Region Size (num elements)');
legend;
title('Feasible Region Size vs. Dimension');
grid on;
hold off;

%% Print summary table
fprintf('\nSummary Table:\n');
for i = 1:length(dims)
    s = all_stats{i};
    fprintf('nx=%d: t_total=%.2fs, t_MRPI=%.2fs, t_feareg=%.2fs, t_sdp=%.2fs, S_hat_opt_vertices=%d\n', ...
        s.nx, s.elapsedTime, s.t_mrpi_shat, s.t_feareg, s.t_sdp, s.S_hat_opt_num_vertices);
end
