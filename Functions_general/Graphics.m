classdef Graphics<handle
    % This is a function collection rather than a class
    % So, All methods are static.
    % plot trajectory and show convex set which are projected from n-dim
    
    methods (Static)
        function show_convex(P, varargin)
            % Allow selection of projection dimensions, default [1 2]
            if nargin > 2 && isnumeric(varargin{1}) && numel(varargin{1}) == 2
                proj_dims = varargin{1};
                varargin = varargin(2:end);
            else
                proj_dims = [1 2];
            end
            P_reduced = Graphics.projectPolytope2Plane(P, proj_dims);
            fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{:})
            hold on;
        end

        function show_trajectory(x_seq, varargin)
            plot(x_seq(1, :), x_seq(2, :), varargin{:})
            hold on;
        end

        function P_projected = projectPolytope2Plane(P, proj_dims)
            % If P is a Zonotope, convert to Polyhedron
            if isa(P, 'Zonotope')
                P = P.toPolyhedron();
            end
            if isprop(P, 'V') && ~isempty(P.V)
                vert = P.V;
            else
                vert = P.vertices();
            end
            x_vert = round(vert(:, proj_dims(1)), 5);
            y_vert = round(vert(:, proj_dims(2)), 5);
            idx = convhull(x_vert, y_vert);
            P_projected = Polyhedron([x_vert(idx), y_vert(idx)]);
        end
    end
end