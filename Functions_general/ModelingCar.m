% ModelingCar.m
%
% Modeling of the ego vehicle (EV) and the leading vehicle (LV).

classdef ModelingCar < handle
    
    properties (SetAccess = public)
        A;
        B;
        Xi_true_EV;
        Xi_true_EV_A;
        Xi_true_EV_b;
        Xi_true_LV;
        Xi_true_LV_A;
        Xi_true_LV_b;
        max_u_LV;
        min_u_LV;
    end
    
    methods (Access = public)
        function obj = ModelingCar(parameters)
            obj.A = parameters.A;
            obj.B = parameters.B;
            obj.Xi_true_EV = parameters.Xi_true_EV;
            obj.Xi_true_EV_A = parameters.Xi_true_EV.A;
            obj.Xi_true_EV_b = parameters.Xi_true_EV.b;
            obj.Xi_true_LV = parameters.Xi_true_LV;
            obj.Xi_true_LV_A = parameters.Xi_true_LV.A;
            obj.Xi_true_LV_b = parameters.Xi_true_LV.b;
            obj.max_u_LV = parameters.max_u_LV;
            obj.min_u_LV = parameters.min_u_LV;
        end

        function [xi_EV_k, x_EV_k_next] = EVModeling(obj, x_EV_k, u_EV_k)
            % system dynamic of the EV
            xi_EV_k = polytope_sample_EV(obj, 1); % random sample from \Xi_EV_true
            x_EV_k_next = obj.A*x_EV_k + obj.B*u_EV_k + xi_EV_k;
        end
        
        function [xi_LV_k, x_LV_k_next] = LVModeling(obj, x_LV_k)
            % system dynamic of the LV
            xi_LV_k = polytope_sample_LV(obj, 1); % random sample from \Xi_LV_true
            u_LV_random = (obj.max_u_LV - obj.min_u_LV)*rand(1) + obj.min_u_LV; % random sample from U_LV_true
            x_LV_k_next = obj.A*x_LV_k + xi_LV_k + obj.B*u_LV_random;
        end
        
        function EV_samples = polytope_sample_EV(obj, N_sam)
            nx = size(obj.Xi_true_EV_A,2);
            V = obj.Xi_true_EV.V;
            v_min = min(V, [], 1);
            v_max = max(V, [], 1);
            EV_samples = zeros(nx, N_sam);
            count = 0;
            max_trials = 1000 * N_sam;
            trials = 0;
            while count < N_sam && trials < max_trials
                x = v_min + (v_max - v_min) .* rand(1, nx);
                if all(obj.Xi_true_EV_A * x' <= obj.Xi_true_EV_b)
                    count = count + 1;
                    EV_samples(:, count) = x';
                end
                trials = trials + 1;
            end
            if count < N_sam
                warning('polytope_sample_EV: Only %d samples generated out of %d requested.', count, N_sam);
                EV_samples = EV_samples(:, 1:count);
            end
        end
        
        function LV_samples = polytope_sample_LV(obj, N_sam)
            nx = size(obj.Xi_true_LV_A,2);
            V = obj.Xi_true_LV.V;
            v_min = min(V, [], 1);
            v_max = max(V, [], 1);
            LV_samples = zeros(nx, N_sam);
            count = 0;
            max_trials = 1000 * N_sam;
            trials = 0;
            while count < N_sam && trials < max_trials
                x = v_min + (v_max - v_min) .* rand(1, nx);
                if all(obj.Xi_true_LV_A * x' <= obj.Xi_true_LV_b)
                    count = count + 1;
                    LV_samples(:, count) = x';
                end
                trials = trials + 1;
            end
            if count < N_sam
                warning('polytope_sample_LV: Only %d samples generated out of %d requested.', count, N_sam);
                LV_samples = LV_samples(:, 1:count);
            end
        end
        function samples = polytope_sample(obj, N_sam, A, b, V)
            % This function is no longer needed for nx > 2, as MPT3 is used
            error('polytope_sample is deprecated: use Polyhedron.uniformSample for nx > 2');
        end
        
    end
end
