function Fs = MRPISet(Ak, W, epsilon)
% Sampling-Based Method for Computing the Robust Positive Invariant Set
[nx, ~] = size(Ak); 
s       = 0; 
alpha   = 1000;
Ms      = 1000;
it      = 0;

%%% CHANGE!
E = eye(nx);
D = [E; -E]; % basis directions
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

Sam_num = 10000;
Sample  = 50*rand(nx,Sam_num) - 25;


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
elseif nx == 3
    [convhull_index,~] = convhull(sample_proj(1,:)', sample_proj(2,:)', sample_proj(3,:)');
else
    error('convhull not implemented for dimensions > 3');
end

MRPI_W             = sample_proj(:,convhull_index);
Fs                 = Polyhedron(MRPI_W');

end
