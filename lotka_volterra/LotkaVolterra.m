classdef LotkaVolterra < handle
    %   This class simulates the classic predator-prey dynamics:
    %       dx/dt = alpha*x - beta*x*y     (prey growth - predation)
    %       dy/dt = delta*x*y - gamma*y    (predator growth - death)
    %
    %   Example usage:
    %       lv = LotkaVolterra();
    %       lv.alpha = 1.0;  % prey growth rate
    %       lv.beta = 0.1;   % predation rate
    %       lv.delta = 0.075; % predator growth rate
    %       lv.gamma = 1.5;  % predator death rate
    %       lv = lv.solve([10; 5]);  % solve with initial conditions
    %       lv.plot_phase_space();
    
    properties
        % default model parameters
        alpha = 1.0  % Prey birth rate
        beta = 0.1  % Predation rate (prey death per predator)
        delta = 0.075  % Predator birth rate (per prey consumed)
        gamma = 1.5   % Predator death rate
        
        % Noise parameter (for stochastic simulation)
        sigma = 0.0     
        
        % Time parameters
        T = 100   % Total simulation time
        dt = 0.01  % Time step
        
        % State variables (set after solving)
        t    % Time vector
        X    % State trajectory [prey, predator] (Nt x 2)
        Nt   % Number of time points
    end
    
    methods
        function obj = LotkaVolterra(opts)
            % Constructor
            % lv = LotkaVolterra() creates with default parameters
            % lv = LotkaVolterra('param', value, ...) sets parameters
            % example lv = LotkaVolterra('T', 50, 'alpha', 2.0)
            arguments
                opts.alpha = 1.0
                opts.beta = 0.1
                opts.delta = 0.075
                opts.gamma = 1.5
                opts.sigma = 0.0
                opts.T = 100
                opts.dt = 0.01
            end
            
            % store given parameter values to obj
            obj.alpha = opts.alpha;
            obj.beta = opts.beta;
            obj.delta = opts.delta;
            obj.gamma = opts.gamma;
            obj.sigma = opts.sigma;
            obj.T = opts.T;
            obj.dt = opts.dt;
        end
        
        function dX_dt = dynamics(obj, X)
            % Compute the deterministic dynamics
            % X = current state of system at ANY point (Need not be initial state)
            % dX_dt = dynamics(obj, X) returns [dx/dt; dy/dt]
            
            x = X(1);  % current prey popln
            y = X(2);  % current predator popln
             
            dx_dt = obj.alpha * x - obj.beta * x * y;
            dy_dt = obj.delta * x * y - obj.gamma * y;
            
            dX_dt = [dx_dt; dy_dt]; % return rate of change
        end
        
        function obj = solve(obj, X0, opts)
            % Find approximate soln to Lotka-Volterra equations using RK4
            % obj = solve(obj, X0) solves from initial condition X0 using RK4
            arguments
                obj   % LV instance
                X0   % initial conditions
                opts.seed = []
            end
            
            if ~isempty(opts.seed)
                rng(opts.seed);
            end
            
            % Initialize
            % notice that if T = 100, dt = 0.01 then Nt = 10001
            obj.Nt = round(obj.T / obj.dt) + 1; % number of points from trajectory to plot
            obj.t = (0:obj.Nt-1)' * obj.dt; % vector of time points
            obj.X = zeros(obj.Nt, 2); % matrix to store prey, predator values at each time point
            obj.X(1, :) = X0(:)'; % store initial condition X0 in first row
            
            % Integrate using Runge-Kutta 4th order
            for i = 1:obj.Nt-1
                Xcurr = obj.X(i, :)';
                
                k1 = obj.dynamics(Xcurr);
                k2 = obj.dynamics(Xcurr + obj.dt/2 * k1);
                k3 = obj.dynamics(Xcurr + obj.dt/2 * k2);
                k4 = obj.dynamics(Xcurr + obj.dt * k3);
                Xnext = Xcurr + obj.dt/6 * (k1 + 2*k2 + 2*k3 + k4);
                
                % Add noise if sigma > 0
                if obj.sigma > 0
                    noise = obj.sigma * sqrt(obj.dt) * randn(2, 1);
                    Xnext = Xnext + noise .* Xnext;  % multiplicative noise
                end
                
                % Ensure non-negative populations
                Xnext = max(Xnext, 0);
                
                obj.X(i+1, :) = Xnext';
            end
        end
        
        function [x_eq, y_eq] = equilibria(obj)
            % Compute the equilibrium points
            %  [x_eq, y_eq] = equilibria(obj) returns equilibrium points
            %
            % There are two equilibria:
            %   1. Trivial: (0, 0)
            %   2. Non-trivial: (gamma/delta, alpha/beta)
            
            % Trivial equilibrium
            x_eq(1) = 0;
            y_eq(1) = 0;
            
            % Non-trivial equilibrium (center of oscillations)
            x_eq(2) = obj.gamma / obj.delta;
            y_eq(2) = obj.alpha / obj.beta;
        end
        
        function fig = plot_time_series(obj)
            %Plot prey and predator populations over time
            
            fig = figure('Name', 'Lotka-Volterra Time Series');
            
            subplot(2,1,1);
            plot(obj.t, obj.X(:,1), 'b-', 'LineWidth', 1.5);
            xlabel('Time');
            ylabel('Prey Population (x)');
            title('Prey Dynamics');
            grid on;
            
            subplot(2,1,2);
            plot(obj.t, obj.X(:,2), 'r-', 'LineWidth', 1.5);
            xlabel('Time');
            ylabel('Predator Population (y)');
            title('Predator Dynamics');
            grid on;
        end
        
        function fig = plot_phase_space(obj, opts)
            % Plot the phase space trajectory
            %   fig = plot_phase_space(obj)
            %   fig = plot_phase_space(obj, 'showEquilibrium', true)
            arguments
                obj
                opts.showEquilibrium = true
                opts.colorByTime = true
            end
            
            fig = figure('Name', 'Lotka-Volterra Phase Space');
            hold on;
            
            if opts.colorByTime
                % Color trajectory by time
                scatter(obj.X(:,1), obj.X(:,2), 10, obj.t, 'filled');
                colorbar;
                colormap(jet);
                cb = colorbar;
                cb.Label.String = 'Time';
            else
                plot(obj.X(:,1), obj.X(:,2), 'b-', 'LineWidth', 1);
            end
            
            % Mark start and end
            plot(obj.X(1,1), obj.X(1,2), 'go', 'MarkerSize', 12, 'LineWidth', 2);
            plot(obj.X(end,1), obj.X(end,2), 'rs', 'MarkerSize', 12, 'LineWidth', 2);
            
            % Show equilibrium
            if opts.showEquilibrium
                [x_eq, y_eq] = obj.equilibria();
                plot(x_eq(2), y_eq(2), 'k+', 'MarkerSize', 15, 'LineWidth', 3);
            end
            
            xlabel('Prey Population (x)');
            ylabel('Predator Population (y)');
            title('Phase Space Trajectory');
            legend('Trajectory', 'Start', 'End', 'Equilibrium', 'Location', 'best');
            grid on;
            axis equal;
            hold off;
        end
        
        % subsampling reduces number of points so faster computation
        function X_sub = subsample(obj, SI)
            %SUBSAMPLE Return subsampled trajectory
            %   X_sub = subsample(obj, SI) subsamples with interval SI
            
            Nds = round(SI / obj.dt);
            idx = 1:Nds:obj.Nt;
            X_sub = obj.X(idx, :);
        end
        
        function t_sub = subsample_time(obj, SI)
            %SUBSAMPLE_TIME Return subsampled time vector
            
            Nds = round(SI / obj.dt);
            idx = 1:Nds:obj.Nt;
            t_sub = obj.t(idx);
        end
    end
end
