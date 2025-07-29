% MRPISet.m
%
% Computes the robust positive invariant set S for a given system.

function Fs = MRPISet(Ak, W, epsilon)
% Sampling-Based Method for Computing the Robust Positive Invariant Set
[nx, ~] = size(Ak); 
s       = 0; 
alpha   = 1000;
Ms      = 1000;
it      = 0;

%%% CHANGE! Add more directions for support function
E = eye(nx);
D = [E; -E]; % axis directions
% Add diagonals
if nx > 2
    D = [D; ones(1,nx)/sqrt(nx); -ones(1,nx)/sqrt(nx)];
end
% Add random directions for more coverage
n_extra = 10;
extra_dirs = randn(n_extra, nx);
extra_dirs = extra_dirs ./ vecnorm(extra_dirs,2,2); % normalize
D = [D; extra_dirs];
D_transpose = D';

while(alpha > epsilon/(epsilon + Ms))
    s     = s + 1;
    alpha = max(W.support(Ak^s*(W.A)')./W.b);
    
    %%% CHANGE!
    mss   = zeros(size(D,1), 1);       % 2*nx rows
    
    %%% CHANGE!
    Ak_powers = cell(s, 1);
    Ak_powers{1} = Ak;
    for i = 2:s
        Ak_powers{i} = Ak_powers{i-1} * Ak;
    end 
    
    dirs = [];
    for i = 1:s
        dirs = [dirs, Ak_powers{i} * D_transpose];  % each direction is a column
    end
    mss_all = W.support(dirs);  % returns column vector
    Ms = max(mss_all);
        
    % for i = 1:s
    %    for j = 1:size(D,1)
    %        dir = (Ak^i) * D(j, :)';  % nx×1
    %        mss(j) = W.support(dir);  % no transpose! MPT3 expects column vector
    %    end
    % end
    
    %for i   = 1:s
    %    %%% CHANGE!
    %    mss = mss + W.support((Ak^i *D')');
    %end
    
    % Ms = max(mss);   
    
    it = it+1;
end

% Sample from the actual set W
Sam_num = 5000;
try
    % If W has a sample method (MPT3 Polyhedron)
    Sample = W.randomPoint(Sam_num)'; % [nx x Sam_num]
catch
    % Fallback: uniform box
    warning('W.randomPoint failed, using uniform box sampling.');
    Sample = 50*rand(nx,Sam_num) - 25;
end


yalmip('clear')
z = sdpvar(nx, s); 
proj_sample  = sdpvar(nx,1);
Input_sample = sdpvar(nx,1);
cns  = [];
item = zeros(nx,1);
for i = 1:s
    item = item + Ak^(i-1)*z(:,i);
    cns  = [cns,W.A*z(:,i) <= W.b];
end
cns  = [cns,proj_sample == item/(1-alpha)];
obj  = norm(proj_sample - Input_sample);
ops  = sdpsettings('relax',0);
MRPI = optimizer(cns,obj,ops, Input_sample, proj_sample);

parfor kkk=1:Sam_num
     [sample_proj(:,kkk), ~] = MRPI(Sample(:,kkk));
end

% [convhull_index,~] = convhull(sample_proj');
if nx == 2
    [convhull_index,~] = convhull(sample_proj');
    MRPI_W = sample_proj(:,convhull_index);
    Fs = Polyhedron(MRPI_W');
elseif nx == 3
    [convhull_index,~] = convhull(sample_proj(1,:)', sample_proj(2,:)', sample_proj(3,:)');
    MRPI_W = sample_proj(:,convhull_index);
    Fs = Polyhedron(MRPI_W');
else
    % For nx > 3, use zonotope approximation with MPT3
    % Center: mean of projected samples
    zono_center = mean(sample_proj, 2);
    % Generators: principal directions (PCA) or random subset of sample differences
    num_gens = min(10, size(sample_proj,2));
    sample_diffs = sample_proj(:,1:num_gens) - zono_center;
    Fs = Zonotope(zono_center, sample_diffs);
end

%{
% --- DEBUG: Print max vertex magnitude of MRPI set after computation ---
% Compute and print the maximum absolute value of any vertex in the MRPI set
if exist('Fs', 'var') && ~isempty(Fs.V)
    max_vertex = max(abs(Fs.V(:)));
    fprintf('DEBUG: Max absolute value in MRPI set vertices: %g\n', max_vertex);
else
    fprintf('DEBUG: MRPI set vertices not available.\n');
end

% --- DEBUG: Check stability of closed-loop matrix Phi if available ---
if exist('Ak', 'var')
    eig_Ak = eig(Ak);
    fprintf('DEBUG: Closed-loop eigenvalues (should be <1 in magnitude for stability):\n');
    disp(eig_Ak);
    if all(abs(eig_Ak) < 1)
        fprintf('DEBUG: Closed-loop system is stable.\n');
    else
        fprintf('DEBUG: WARNING! Closed-loop system is NOT stable.\n');
    end
end

%}
end
