%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Results = PSC(grn, meshxt, ic, psc)

Results.t    = [];
Results.u    = [];
Results.cost = [];
Results.PX   = {};

%% ===================== CONTROL DIMENSION =====================

% Number of inductors
psc.n_I = numel(psc.alpha);

% Check inputs
checkPSCinputs(grn, meshxt, ic, psc);

% Infer inductor → gene mapping (automatic)
psc.inducer2gene = inferInducerTargets(grn, meshxt, psc);

% Size of PSC search space
psc.size_psc = 2^psc.n_I;

% All possible inducer configurations
psc.S = dec2bin(0:psc.size_psc-1) - '0';

%% ===================== KAPPA COMPUTATION =====================

psc.kappa = kappa_value(grn, meshxt, psc);

% Control input matrix
psc.u = psc.kappa .* psc.S;

%% ===================== TIME GRID =====================

tloop = linspace(meshxt.Time_mesh(1), ...
                 meshxt.Time_mesh(2), ...
                 meshxt.Time_mesh(3)+1);

cloop = tloop(1:meshxt.w:end);

Results.t = cloop(1);

%% ===================== INITIAL CONDITION =====================

solution_i = cell(1, psc.size_psc);
cost_i = zeros(1, psc.size_psc);

PX0 = initial_condition(grn, meshxt, ic);

Results.PX   = {PX0};
Results.cost = psc.cost_f(PX0);

%% ===================== MAIN PSC LOOP =====================

for t = 1:length(cloop)-1
    
    for r = 1:psc.size_psc
        
        % Solve PIDE
        solution_i{r} = PIDE(grn, meshxt, psc, t, r, PX0);
        PX         = solution_i{r}.PTX;

        cost_i(r) = psc.cost_f(PX);
    end

    % Select optimal control
    if psc.opt == 1
        [cost, ind] = max(cost_i);
    else
        [cost, ind] = min(cost_i);
    end

    u_t = psc.u(ind, :);

    % Propagate with selected control
    solution   = solution_i{ind};
    PX         = solution.PTX;
    PX0        = gather(PX);

    Results.t    = [Results.t; cloop(t+1)];
    Results.u    = [Results.u; u_t];
    Results.cost = [Results.cost; cost];
end

Results.PX = [Results.PX; {PX}];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ind2gene = inferInducerTargets(grn, meshxt, psc)

n_I    = psc.n_I;
n_gene = grn.n_gene;

% Max protein values
max_x = cell(1, n_gene);
for g = 1:n_gene
    max_x{g} = meshxt.Prot_mesh(g,2);
end

ind2gene = zeros(1, n_I);

tol = 1e-10;

for i = 1:n_I
    
    I0 = zeros(1, n_I);
    I1 = zeros(1, n_I);
    I1(i) = 1;   % arbitrary non-zero activation
    
    for g = 1:n_gene
        
        try
            r0 = grn.rho{g}(max_x, I0);
            r1 = grn.rho{g}(max_x, I1);
        catch
            continue
        end
        
        if abs(r1 - r0) > tol
            ind2gene(i) = g;
            break
        end
    end
    
    if ind2gene(i) == 0
        error('Inducer %d does not affect any rho{g}. Model ill-defined.', i);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kappa = kappa_value(grn, meshxt, psc)

n_I    = psc.n_I;
n_gene = grn.n_gene;

kappa = NaN(1, n_I);

% Max protein values
max_x = cell(1, n_gene);
for g = 1:n_gene
    max_x{g} = meshxt.Prot_mesh(g,2);
end

% Inducer search grid
I_vals = logspace(-4, 10, 1e4);

for i = 1:n_I
    
    alpha_i = psc.alpha(i);
    g_i     = psc.inducer2gene(i);  
    
    % rho_g with only inducer i active
    rho_i = @(u) grn.rho{g_i}( ...
        max_x, ...
        build_I_vector(u, i, n_I) );
    
    rho_vals = arrayfun(rho_i, I_vals);
    target   = 1 - alpha_i;

    diffs = rho_vals - target;
    idx   = find(diffs(1:end-1).*diffs(2:end) < 0, 1);

    if isempty(idx)
        warning('No kappa found for inducer %d (gene %d).', i, g_i);
        continue;
    end

    kappa(i) = fzero(@(u) rho_i(u) - target, ...
                     [I_vals(idx), I_vals(idx+1)]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function I = build_I_vector(u, i_active, n_I)

I = zeros(1, n_I);
I(i_active) = u;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PX0 = initial_condition(grn, meshxt, ic)

% Spatial discretization
if gpuDeviceCount > 0 
    meshxt.Prot_mesh_gpu = gpuArray(meshxt.Prot_mesh);
    
    iN = cell(grn.n_gene,1); 
    x = cell(grn.n_gene,1);
    Xgrid = cell(grn.n_gene,1);
    for i = 1:grn.n_gene
        iN{i} = meshxt.Prot_mesh_gpu(i,3) + 1;
        x{i} = linspace(meshxt.Prot_mesh_gpu(i,1), meshxt.Prot_mesh_gpu(i,2), iN{i}); 
    end
    [Xgrid{1:grn.n_gene}] = ndgrid(x{:});

else 

    iN=cell(grn.n_gene,1);
    x=cell(grn.n_gene,1);
    for i=1:grn.n_gene
        iN{i} = meshxt.Prot_mesh(i,3) + 1;
        x{i} = linspace(meshxt.Prot_mesh(i,1),meshxt.Prot_mesh(i,2), iN{i});
    end
    Xgrid=cell(grn.n_gene,1);
    [Xgrid{1:grn.n_gene}] = ndgrid(x{1:grn.n_gene});
end

if isfield(ic,'PX0')
    % OPTION 1: load distribution
    PX0 = ic.PX0;
else
    % OPTION 0: create Gaussian density function using mean and variance
    gausker = 0;
    for i = 1:grn.n_gene
        gausker = gausker + ((Xgrid{i} - ic.mean(i)).^2) / ic.var(i);
    end
    PX0_aux = exp(-gausker/2);
    
    % Norm of the initial condition
    auxnor0=PX0_aux;
    for i=1:grn.n_gene
        auxnor0 = trapz(x{i},auxnor0);
    end
    
    % Normalization of the initial condition to be a density function
    PX0=PX0_aux/auxnor0;

end
    if gpuDeviceCount > 0 
        PX0 = gpuArray(PX0);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkPSCinputs(grn, meshxt, ic, psc)

% CHECKPSCINPUTS - Checks consistency and validity of PSC input parameters
%
%   This function throws errors or warnings if any of the PSC inputs
%   are inconsistent. Recommended to run at the beginning of PSC or main.

% ====== 1. Basic GRN ======
assert(isstruct(grn), 'GRN must be a struct.');
required_fields = {'n_gene','k_m','k_x','g_m','g_x','eps','rho','F'};
for f = required_fields
    if ~isfield(grn,f{1})
        error('GRN missing field: %s', f{1});
    end
end

% Check dimensions
n = grn.n_gene;
vec_fields = {'k_m','k_x','g_m','g_x','eps'};
for f = vec_fields
    if length(grn.(f{1})) ~= n
        error('GRN.%s must have length equal to n_gene (%d)', f{1}, n);
    end
end

% Check function cells
if ~iscell(grn.rho) || length(grn.rho) ~= n
    error('GRN.rho must be a cell array of length n_gene.');
end
if ~iscell(grn.F) || length(grn.F) ~= psc.n_I
    error('GRN.F must be a cell array of length n_gene.');
end

% ====== 2. meshxt ======
assert(isstruct(meshxt),'meshxt must be a struct.');
required_meshxt_fields = {'Prot_mesh','Time_mesh','w'};
for f = required_meshxt_fields
    if ~isfield(meshxt,f{1})
        error('meshxt missing field: %s', f{1});
    end
end

% Check Prot_mesh dimensions
if size(meshxt.Prot_mesh,1) ~= n
    error('meshxt.Prot_mesh rows (%d) must match n_gene (%d).', size(meshxt.Prot_mesh,1), n);
end

% Check Time_mesh
if length(meshxt.Time_mesh) ~= 3
    error('meshxt.Time_mesh must have 3 elements: [t0 dt t_final].');
end

% Check control window
if meshxt.w < 1 || mod(meshxt.w,1) ~= 0
    warning('meshxt.w should be a positive integer. Adjusted to 1.');
    meshxt.w = 1;
end

% Check if w exceeds the simulation horizon
n_total_points = meshxt.Time_mesh(3);
if meshxt.w > n_total_points
    error('Inconsistent Window Size: meshxt.w (%d) cannot be larger than the total time points (%d).', ...
        meshxt.w, n_total_points);
end

% ====== 3. Initial Condition ======
assert(isstruct(ic), 'ic must be a struct.');
assert(isfield(ic,'mode'), 'ic must define ic.mode.');
assert(ismember(ic.mode,[0 1]), 'ic.mode must be 0 (Gaussian) or 1 (PX0).');

switch ic.mode

    case 0
        % Gaussian IC
        assert(isfield(ic,'mean') && isfield(ic,'var'), ...
            'Gaussian IC requires ic.mean and ic.var.');

        assert(length(ic.mean) == n && length(ic.var) == n, ...
            'ic.mean and ic.var must have length equal to n_gene.');

    case 1
        % External distribution
        assert(isfield(ic,'PX0') && ~isempty(ic.PX0), ...
            'PX0 mode requires a valid ic.PX0.');

end

% ====== 4. PSC ======
assert(isstruct(psc),'psc must be a struct.');
if  ~isfield(psc,'cost_f') || ~isfield(psc,'opt')
    error('psc must have cost_f, and opt fields.');
end

% Consistency checks
if ~iscell(grn.F) || numel(grn.F) ~= psc.n_I
error('PSC:InducerMismatch', ...
'Number of inducer functions grn.F must match length of psc.alpha.');
end

% Cost function
if ~isa(psc.cost_f,'function_handle')
    error('psc.cost_f must be a function handle.');
end
try
    test_val = psc.cost_f(ones([meshxt.Prot_mesh(:,3)' + 1, 1]));
    if ~isscalar(test_val)
        warning('psc.cost_f should return a scalar value. Returned value might be ignored.');
    end
catch ME
    error('psc.cost_f failed during test: %s', ME.message);
end


% ====== 5. General Checks ======
% Negative values
if any(grn.k_m<0) || any(grn.k_x<0) || any(grn.g_m<0) || any(grn.g_x<0)
    warning('Some GRN rates are negative. Please check units and values.');
end

%disp('PSC input check passed successfully.');
end

