%% Temporal Mapper for Lotka-Volterra ODE
% Temporal mapper construction using the Lotka-Volterra predator-prey model 
% A simpler system to understand temporal method.
%
% The temporal mapper workflow:
% 1. Simulate/obtain time series data
% 2. Construct a time-ordered k-nearest neighbor graph (tknndigraph.m)
% 3. Filter the graph to create a "shape graph" (filtergraph.m)
% 4. Analyze: geodesic distances, recurrence plots (TCM)
%
% Based on tmapper repository (Mengsen Zhang): https://github.com/braindynamicslab/tmapper

clear all
close all
clc

% Add the main tmapper code to the path (relative to this script's location)
addpath('./code/');


%% Simulating Lotka-Volterra System (LotkaVolterra.m)



% Create model with parameters that produce oscillations
lv = LotkaVolterra();
lv.alpha = 1;   % prey birth rate 1.0
lv.beta = 0.1;      % predation rate 0.1
lv.delta = 0.075;    % predator birth rate from consuming prey 0.075
lv.gamma = 1;     % predator death rate 1.5
lv.T = 200;       % simulation time 
lv.dt = 0.01;    % integration time step

% Initial conditions (away from equilibrium to see oscillations)
X0 = [10; 5];

% Solve the system
fprintf('Solve ODE using RK4 given the initial condition: prey=%.1f, predator=%.1f\n', X0(1), X0(2));
lv = lv.solve(X0);

% Display equilibrium point
[x_eq, y_eq] = lv.equilibria();
fprintf('Non-trivial equilibrium (center of oscillations): prey=%.2f, predator=%.2f\n\n', x_eq(2), y_eq(2));

% Plot the raw dynamics
figure('Name', 'Lotka-Volterra Dynamics', 'Position', [100 100 600 900]);

subplot(3,1,1);
plot(lv.t, lv.X(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(lv.t, lv.X(:,2), 'r-', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Population');
title('Predatory-Prey time series data');
legend('Prey', 'Predator');
grid on;

subplot(3,1,2);
plot(lv.X(:,1), lv.X(:,2), 'k-', 'LineWidth', 0.5);
hold on;
scatter(lv.X(:,1), lv.X(:,2), 5, lv.t, 'filled');
plot(x_eq(2), y_eq(2), 'g+', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Prey'); ylabel('Predator');
title('Phase Space (colored by time)');
colorbar; colormap(jet);
grid on;

subplot(3,1,3);
% Show the cyclic nature - prey vs time colored by predator
scatter(lv.t, lv.X(:,1), 5, lv.X(:,2), 'filled');
xlabel('Time'); ylabel('Prey Population');
title('Prey (colored by predator)');
colorbar; colormap(jet);
grid on;

exportgraphics(gcf, 'images/lv_dynamics.pdf', 'ContentType', 'vector');

%% Subsampling Data 

% Subsample interval (like sampling interval in neural recordings)
% if SI = 0.5 sec (i.e. desired sampling interval) while original time step (lv.dt) = 0.01
% Nds = round(0.5/0.01) = 50 i.e. take every 50 th point
% spidx = 
SI = 0.5;  % seconds
Nds = round(SI / lv.dt); % num. of downsample steps
spidx = 1:Nds:lv.Nt;    % subsample indices
N_t = length(spidx); % num. of subsampled points


% This is like the temporal resolution of fMRI or other recordings.
fprintf('Original time points: %d\n', lv.Nt);
fprintf('Subsampled time points: %d (every %.2f seconds)\n', N_t, SI);


% Extract subsampled data
x_sub = lv.X(spidx, :);
t_sub = lv.t(spidx);
tidx = (1:N_t)';  % Time indices for temporal neighbors


%% Construct Temporal k-NN Graph (../code/tknndigraph.m)


% Parameters for temporal mapper
k = 3;   % number of nearest neighbors in state space
d = 2;   % distance threshold for filtration

% The key insight of temporal mapper:
% - We find k-nearest neighbors in state space i.e. prey-predator coordinates
% - But we ALSO connect consecutive time points (temporal neighbors)
% - This captures both the geometry of the attractor AND temporal ordering

fprintf('Building temporal k-NN graph...\n');
fprintf('This graph connects:\n');
fprintf('  1. Each point to its %d nearest neighbors in state space\n', k);
fprintf('  2. Each point to its temporal neighbors (t -> t+1)\n\n');

g = tknndigraph(x_sub, k, tidx, 'reciprocal', true, 'timeExcludeSpace', true);

% Visualize the temporal k-NN graph
figure('Name', 'Temporal k-NN Graph', 'Position', [100 100 800 600]);
h = plot(g, 'Layout', 'force');
h.NodeCData = tidx;     % color by time index
h.MarkerSize = 4;
colormap(jet);
colorbar;
title(sprintf('Temporal k-NN Graph (k=%d)\n%d subsampled points, %d nodes, %d edges', k, N_t, numnodes(g), numedges(g)));
xlabel('Colored by time');

% Save as PDF
exportgraphics(gcf, 'images/tknn_graph.pdf', 'ContentType', 'vector');


%% Filtering to Create Shape Graph (../code/filtergraph.m)

% The filtration step:
% - Nodes within distance d are collapsed together
% - This simplifies the graph while preserving topological structure
% - Similar to how Mapper creates a "shape graph"

fprintf('Filtering graph with threshold d=%d...\n', d);
fprintf('Nodes that are mutual neighbors i.e. both nodes are within <%d steps are merged into one new node.\n\n', d);

[g_simp, members, nsize] = filtergraph(g, d, 'reciprocal', true);

fprintf('Simplified graph: %d nodes, %d edges\n', numnodes(g_simp), numedges(g_simp));
fprintf('Node sizes (number of time points per node):\n');
disp(nsize');


% Create a dedicated figure for the simplified shape graph
figure('Name', 'Simplified Shape Graph (g_simp)', 'Position', [100 100 1200 600]);

% Left subplot: Force-directed layout with node sizes
subplot(1,2,1);
h1 = plot(g_simp, 'Layout', 'force');
h1.NodeCData = 1:numnodes(g_simp);
h1.MarkerSize = sqrt(nsize) * 4;
h1.LineWidth = 1.5;
colormap(gca, jet);
cbar1 = colorbar;
cbar1.Label.String = 'Node Index';
title(sprintf('Shape Graph - Force Layout\n%d nodes, %d edges', numnodes(g_simp), numedges(g_simp)));
xlabel('X'); ylabel('Y');
grid on;

% Right subplot: Circular layout for better visibility of cyclic structure
subplot(1,2,2);
h2 = plot(g_simp, 'Layout', 'circle');
h2.NodeCData = 1:numnodes(g_simp);
h2.MarkerSize = sqrt(nsize) * 4;
h2.LineWidth = 1.5;
colormap(gca, jet);
cbar2 = colorbar;
cbar2.Label.String = 'Node Index';
title(sprintf('Shape Graph - Circular Layout\n(Easier to see cyclic structure)'));
xlabel('X'); ylabel('Y');
grid on;

% Add node information text
fprintf('\nShape Graph Details:\n');
fprintf('Number of nodes: %d\n', numnodes(g_simp));
fprintf('Number of edges: %d\n', numedges(g_simp));
fprintf('Node sizes (time points per node):\n');
for i = 1:numnodes(g_simp)
    fprintf('  Node %d: %d time points\n', i, nsize(i));
end

exportgraphics(gcf, 'images/shape_graph.pdf', 'ContentType', 'vector');

%% Compute Geodesic Distances (../code/normgeo.m)

% Geodesic distance = shortest path length on the graph
% This captures how "far apart" different states are in terms of the dynamics

geod = distances(g_simp, 'Method', 'unweighted');
[geod_n, mNode] = normgeo(geod, nsize); %[distance matrix, new node sizes]

fprintf('Geodesic distance matrix size: %d x %d\n', size(geod));
fprintf('Maximum geodesic distance: %d\n', max(geod(~isinf(geod))));
fprintf('\nGeodesic distances tell us the "dynamical distance" between states.\n');
fprintf('States that are close in state space but far in the dynamics\n');
fprintf('(e.g., different phases of the cycle) will have large geodesic distance.\n\n');



%% Compute Temporal Connectivity Matrix and Recurrence Plot (../code/TCMdistance.m ; ../code/normtcm.m) 

% The TCM shows the shortest path length between all pairs of time points
% This is a form of "recurrence plot" based on graph geodesics

tcm = TCMdistance(g_simp, members); % compute temporal connectivity matrix (tcm)
tcm_n = normtcm(tcm, 'normtype', 'max'); % normalize tcm using 'norm' or 'max' options for 'normtype'

figure('Name', 'Temporal Connectivity Matrix', 'Position', [100 100 800 600]);
imagesc(t_sub, t_sub, tcm_n);
colorbar;
xlabel('Time'); ylabel('Time');
title('Temporal Connectivity Matrix (Recurrence Plot)');
axis square;

exportgraphics(gcf, 'images/tcm.pdf', 'ContentType', 'vector');

fprintf('TCM size: %d x %d (one entry per time point pair)\n', size(tcm));
fprintf('\nThe TCM reveals temporal structure:\n');
fprintf('- Diagonal = 0 (same time point)\n');
fprintf('- Off-diagonal shows how "far" two time points are in the shape graph\n');
fprintf('- Periodic structure appears as repeated patterns\n\n');

%% Visualize Results 

% Create comprehensive visualization
figure('Name', 'Temporal Mapper Results', 'Position', [50 50 1400 800]);

% 7a. Phase space with node coloring
subplot(2,3,1);
% Assign each time point to its node
node_assignment = zeros(N_t, 1);
for i = 1:length(members)
    node_assignment(members{i}) = i;
end
scatter(x_sub(:,1), x_sub(:,2), 30, node_assignment, 'filled');
xlabel('Prey'); ylabel('Predator');
title(sprintf('Phase Space (colored by mapper node)\n%d nodes', length(members)));
colorbar; colormap(gca, jet);
grid on;

% 7b. Time series colored by node
subplot(2,3,2);
scatter(t_sub, x_sub(:,1), 30, node_assignment, 'filled');
xlabel('Time'); ylabel('Prey Population');
title('Prey Time Series (colored by node)');
colorbar; colormap(gca, jet);
grid on;

% 7c. Geodesic distance matrix
subplot(2,3,3);
imagesc(geod);
colorbar;
xlabel('Node'); ylabel('Node');
title('Geodesic Distance Matrix');
axis square;

% 7d. Node sizes (histogram)
subplot(2,3,4);
bar(nsize);
xlabel('Node Index'); ylabel('Number of Time Points');
title('Node Sizes (Points per Node)');
grid on;

% 7e. Temporal Connectivity Matrix (Recurrence Plot)
subplot(2,3,5);
imagesc(t_sub, t_sub, tcm_n);
colorbar;
xlabel('Time'); ylabel('Time');
title('Temporal Connectivity Matrix (Recurrence Plot)');
axis square;

% 7f. Shape graph visualization
subplot(2,3,6);
h = plot(g_simp, 'Layout', 'force');
h.NodeCData = 1:numnodes(g_simp);
h.MarkerSize = sqrt(nsize) * 3;
colormap(gca, jet);
title(sprintf('Shape Graph\n(%d nodes, %d edges)', numnodes(g_simp), numedges(g_simp)));


%% Check different parameter values

% Show how k and d affect the result
figure('Name', 'Parameter Effects', 'Position', [100 100 1000 1000]);

k_values = [2, 3, 4, 5];
for i = 1:4
    k_test = k_values(i);
    g_test = tknndigraph(x_sub, k_test, tidx, 'reciprocal', true, 'timeExcludeSpace', true);
    [g_simp_test, members_test, nsize_test] = filtergraph(g_test, d, 'reciprocal', true);
    
    subplot(2,2,i);
    h = plot(g_simp_test, 'Layout', 'force');
    h.MarkerSize = sqrt(nsize_test) * 2;
    title(sprintf('k=%d: %d nodes', k_test, numnodes(g_simp_test)));
end
sgtitle('Effect of k (number of nearest neighbors) on Shape Graph');

fprintf('Observations:\n');
fprintf('- Smaller k: More nodes, finer resolution of dynamics\n');
fprintf('- Larger k: Fewer nodes, coarser but more robust structure\n');
fprintf('- The cyclic topology should be preserved across reasonable k values\n\n');

%% Save the workspace

outfile = 'lotka_volterra_results.mat';
save(outfile, 'lv', 'x_sub', 't_sub', 'g', 'g_simp', 'members', 'nsize', ...
     'geod', 'geod_n', 'tcm', 'tcm_n', 'k', 'd', 'SI');


