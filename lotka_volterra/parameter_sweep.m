%% Parameter Sweep: Explore Temporal Mapper with Lotka-Volterra
% This script explores how temporal mapper results change with different
% dynamical regimes and mapper parameters.
%
% Useful for understanding:
% - How k (nearest neighbors) affects resolution
% - How d (filtration threshold) affects graph structure  
% - How dynamical parameters change the mapper output

clear all
close all
clc

addpath('./code');

%% ===== EXPERIMENT 1: Effect of initial conditions =====
% Lotka-Volterra conserves a quantity, so different initial conditions
% give different "energy levels" (orbit sizes)

fprintf('=== Experiment 1: Different Initial Conditions ===\n\n');

figure('Name', 'Effect of Initial Conditions', 'Position', [50 50 1500 500]);

initial_conditions = {[5; 5], [10; 5], [20; 10]};
titles = {'Small orbit', 'Medium orbit', 'Large orbit'};

for exp = 1:3
    lv = LotkaVolterra('T', 100, 'dt', 0.01);
    lv = lv.solve(initial_conditions{exp});
    
    % Subsample and build mapper
    SI = 0.5;
    Nds = round(SI / lv.dt);
    spidx = 1:Nds:lv.Nt;
    x_sub = lv.X(spidx, :);
    tidx = (1:length(spidx))';
    
    k = 8; d = 4;
    g = tknndigraph(x_sub, k, tidx, 'reciprocal', true, 'timeExcludeSpace', true);
    [g_simp, members, nsize] = filtergraph(g, d, 'reciprocal', true);
    
    subplot(1,3,exp);
    scatter(x_sub(:,1), x_sub(:,2), 20, tidx, 'filled');
    hold on;
    [x_eq, y_eq] = lv.equilibria();
    plot(x_eq(2), y_eq(2), 'k+', 'MarkerSize', 15, 'LineWidth', 3);
    xlabel('Prey'); ylabel('Predator');
    title(sprintf('%s\nIC: [%.0f, %.0f], Nodes: %d', ...
          titles{exp}, initial_conditions{exp}(1), initial_conditions{exp}(2), ...
          numnodes(g_simp)));
    colorbar; colormap(jet);
    grid on; axis equal;
end

sgtitle('Effect of Initial Conditions on Limit Cycle Size');

%% ===== EXPERIMENT 2: Effect of k (nearest neighbors) =====
fprintf('\n=== Experiment 2: Effect of k (Nearest Neighbors) ===\n\n');

% Setup base simulation
lv = LotkaVolterra('T', 150, 'dt', 0.01);
lv.alpha = 1.0;   % prey birth rate
lv.beta = 0.5;    % predation rate
lv.gamma = 0.8;   % predator death rate
lv.delta = 0.4;   % predator reproduction from prey
lv = lv.solve([10; 5]);

SI = 0.5;
Nds = round(SI / lv.dt);
spidx = 1:Nds:lv.Nt;
x_sub = lv.X(spidx, :);
tidx = (1:length(spidx))';

figure('Name', 'Effect of k', 'Position', [50 100 1500 800]);

k_values = [2, 3, 4, 5, 6, 7];
d = 3;

for i = 1:6
    k = k_values(i);
    g = tknndigraph(x_sub, k, tidx, 'reciprocal', true, 'timeExcludeSpace', true);
    [g_simp, members, nsize] = filtergraph(g, d, 'reciprocal', true);
    
    % TCM
    tcm = TCMdistance(g_simp, members);
    tcm_n = normtcm(tcm, 'normtype', 'max');
    
    % Plot shape graph
    subplot(2,6,i);
    h = plot(g_simp, 'Layout', 'force');
    h.MarkerSize = max(3, sqrt(nsize) * 2);
    h.NodeCData = 1:numnodes(g_simp);
    colormap(gca, jet);
    title(sprintf('k=%d\n%d nodes', k, numnodes(g_simp)));
    
    % Plot TCM
    subplot(2,6,i+6);
    imagesc(tcm_n);
    axis square;
    title('TCM');
    colorbar;
end

sgtitle(sprintf('Effect of k (nearest neighbors), d=%d fixed', d));

fprintf('Key insight: Larger k -> Fewer nodes (more points are neighbors)\n');
fprintf('             Smaller k -> More nodes (finer resolution)\n\n');

%% ===== EXPERIMENT 3: Effect of d (filtration threshold) =====
fprintf('=== Experiment 3: Effect of d (Filtration Threshold) ===\n\n');

figure('Name', 'Effect of d', 'Position', [50 100 1500 800]);

k = 8;
d_values = [2, 3, 4, 6, 8, 12];

for i = 1:6
    d = d_values(i);
    g = tknndigraph(x_sub, k, tidx, 'reciprocal', true, 'timeExcludeSpace', true);
    [g_simp, members, nsize] = filtergraph(g, d, 'reciprocal', true);
    
    % TCM
    tcm = TCMdistance(g_simp, members);
    tcm_n = normtcm(tcm, 'normtype', 'max');
    
    % Plot shape graph
    subplot(2,6,i);
    h = plot(g_simp, 'Layout', 'force');
    h.MarkerSize = max(3, sqrt(nsize) * 2);
    h.NodeCData = 1:numnodes(g_simp);
    colormap(gca, jet);
    title(sprintf('d=%d\n%d nodes', d, numnodes(g_simp)));
    
    % Plot TCM
    subplot(2,6,i+6);
    imagesc(tcm_n);
    axis square;
    title('TCM');
    colorbar;
end

sgtitle(sprintf('Effect of d (filtration threshold), k=%d fixed', k));

fprintf('Key insight: Larger d -> Fewer nodes (more merging)\n');
fprintf('             Smaller d -> More nodes (less merging)\n\n');

%% ===== EXPERIMENT 4: Effect of noise =====
fprintf('=== Experiment 4: Effect of Noise ===\n\n');

figure('Name', 'Effect of Noise', 'Position', [50 100 1500 500]);

sigma_values = [0, 0.01, 0.05];
titles = {'No noise', '\sigma=0.01', '\sigma=0.05'};

for exp = 1:3
    lv = LotkaVolterra('T', 100, 'dt', 0.01, 'sigma', sigma_values(exp));
    lv = lv.solve([10; 5], 'seed', 42);
    
    % Subsample and build mapper
    SI = 0.5;
    Nds = round(SI / lv.dt);
    spidx = 1:Nds:lv.Nt;
    x_sub = lv.X(spidx, :);
    tidx = (1:length(spidx))';
    
    k = 8; d = 4;
    g = tknndigraph(x_sub, k, tidx, 'reciprocal', true, 'timeExcludeSpace', true);
    [g_simp, members, nsize] = filtergraph(g, d, 'reciprocal', true);
    
    subplot(1,3,exp);
    scatter(x_sub(:,1), x_sub(:,2), 20, tidx, 'filled');
    xlabel('Prey'); ylabel('Predator');
    title(sprintf('%s\nNodes: %d', titles{exp}, numnodes(g_simp)));
    colorbar; colormap(jet);
    grid on;
end

sgtitle('Effect of Noise on Trajectories');

fprintf('Key insight: Noise can cause trajectory to drift between orbits,\n');
fprintf('             potentially revealing different basins of attraction.\n\n');

%% ===== EXPERIMENT 5: Multiple periods comparison =====
fprintf('=== Experiment 5: Effect of Simulation Length ===\n\n');

figure('Name', 'Effect of Simulation Length', 'Position', [50 100 1500 800]);

T_values = [50, 100, 200, 400];
k = 8; d = 4;

for i = 1:4
    lv = LotkaVolterra('T', T_values(i), 'dt', 0.01);
    lv = lv.solve([10; 5]);
    
    SI = 0.5;
    Nds = round(SI / lv.dt);
    spidx = 1:Nds:lv.Nt;
    x_sub = lv.X(spidx, :);
    tidx = (1:length(spidx))';
    
    g = tknndigraph(x_sub, k, tidx, 'reciprocal', true, 'timeExcludeSpace', true);
    [g_simp, members, nsize] = filtergraph(g, d, 'reciprocal', true);
    
    % TCM
    tcm = TCMdistance(g_simp, members);
    tcm_n = normtcm(tcm, 'normtype', 'max');
    
    % Plot shape graph
    subplot(2,4,i);
    h = plot(g_simp, 'Layout', 'force');
    h.MarkerSize = max(3, sqrt(nsize) * 2);
    h.NodeCData = 1:numnodes(g_simp);
    colormap(gca, jet);
    title(sprintf('T=%d, N=%d\n%d nodes', T_values(i), length(spidx), numnodes(g_simp)));
    
    % Plot TCM
    subplot(2,4,i+4);
    imagesc(tcm_n);
    axis square;
    title('TCM');
    colorbar;
end

sgtitle('Effect of Simulation Length (more periods)');

fprintf('Key insight: Longer simulations reveal periodicity more clearly in TCM.\n');
fprintf('             Shape graph stabilizes once cycle is well-sampled.\n\n');

%% ===== Summary =====
fprintf('=== PARAMETER SWEEP SUMMARY ===\n\n');
fprintf('For temporal mapper analysis:\n\n');
fprintf('1. k (nearest neighbors):\n');
fprintf('   - Start with k ~ 5-10%% of N_timepoints\n');
fprintf('   - Too small: noisy, too many nodes\n');
fprintf('   - Too large: loses resolution\n\n');
fprintf('2. d (filtration threshold):\n');
fprintf('   - Start with d ~ 3-5 (path length units)\n');
fprintf('   - Too small: many disconnected components\n');
fprintf('   - Too large: everything collapses to few nodes\n\n');
fprintf('3. Sampling interval (SI):\n');
fprintf('   - Should capture dynamics adequately\n');
fprintf('   - ~10-20 samples per oscillation period is good\n\n');
fprintf('4. Simulation length:\n');
fprintf('   - At least 2-3 complete cycles for periodic systems\n');
fprintf('   - More cycles improve TCM periodicity detection\n');
