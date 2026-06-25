%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solution = PIDE(grn, meshxt, psc, t, r, PX0)

% Dimensionless parameters
b = cell(grn.n_gene,1); 
for i = 1:grn.n_gene
    b{i} = grn.k_x(i)/grn.g_m(i); 
end

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

% Time definition
deltat = (meshxt.Time_mesh(2) - meshxt.Time_mesh(1))/meshxt.Time_mesh(3);
t0     = meshxt.Time_mesh(1) + (t-1)*meshxt.w*deltat;
tmax   = meshxt.Time_mesh(1) + t*meshxt.w*deltat;
nt     = meshxt.w+1;
tl     = linspace(t0, tmax, nt);


% Computation of characteristics curves
xbar = cell(grn.n_gene,1);
xbarlim=cell(grn.n_gene,1);
for i = 1:grn.n_gene
    xbar{i} = x{i}*exp(deltat*grn.g_x(i)); 
    xbarlim{i} = find(xbar{i}>=x{i}(end));
end

% Caracteristics grid
Xbargrid=cell(grn.n_gene,1);
[Xbargrid{1:grn.n_gene}] = ndgrid(xbar{1:grn.n_gene});

% Initialization 
PX = PX0;

% Time Independent functions
e_x = cell(grn.n_gene,1);
e_lx = cell(grn.n_gene,1);
for i = 1:grn.n_gene
    e_x{i} = exp(Xgrid{i}/b{i});
    e_lx{i} = exp(-Xgrid{i}/b{i});
    e_x{i}(isinf(e_x{i})) = realmax;
    e_lx{i}(isinf(e_lx{i})) = realmax;
end

% Input function c
cx = cell(1,grn.n_gene);

u = psc.u(r,:);

if gpuDeviceCount > 0
    Xgrid_gpu = cell(grn.n_gene,1);
    for i = 1:grn.n_gene
        Xgrid_gpu{i} = gpuArray(Xgrid{i}); 
    end

    for i = 1:grn.n_gene
    cx{i} = grn.c{i}(Xgrid_gpu,u);
    end
else
    for i = 1:grn.n_gene
    cx{i} = grn.c{i}(Xgrid,u);
    end
end


% Other time independent functions
sumkmcx = 0; sumprotdeg = 0;
for i = 1:grn.n_gene 
    sumkmcx = sumkmcx + grn.k_m(i)*cx{i}; 
    sumprotdeg = sumprotdeg + grn.g_x(i);
end
expl_den = 1+(sumkmcx-sumprotdeg)*deltat;

% Saving
PX_sol      = [];
TT          = [];

for j = 2:nt

    % PX_bar construction using interpolation 
    if grn.n_gene==1
        PX_bar = interp1(x{1}, PX, Xbargrid{1}); % x{1} in case of regular mesh (linspace), in case of irregular: Xgrid{1}
        PX_bar(isnan(PX_bar)) = 0;
    else
        PX_bar = interpn(x{:}, PX, Xbargrid{:}); % x{:} in case of regular mesh (linspace), in case of irregular: Xgrid{:}
        PX_bar(isnan(PX_bar)) = 0;
    end

    % Integral term computation by numerical integration
    Lix = 0;
    for i = 1:grn.n_gene
        Lix = Lix + grn.k_m(i)/b{i} .* e_lx{i} .* cumtrapz(x{i}, e_x{i} .* cx{i} .* PX, i); 
    end

    % Explicit method    
    PX = (PX_bar+deltat*Lix)./expl_den;

    % Zero boundary condition
    CFaux=cell(grn.n_gene,1);
    for i=1:grn.n_gene
        CFaux{i}=':';
    end
    for i = 1:grn.n_gene
        CF = CFaux; 
        CF{i} = iN{i}; 
        PX(CF{:}) = zeros(size(PX(CF{:}))); 
    end

    % Normalization: int_xmin^xmax(PX)dx=1
    auxnorpx = PX; 
    for i = 1:grn.n_gene
        auxnorpx = trapz(x{i}, auxnorpx); 
    end

    PX = PX/auxnorpx;

    % Saving the solution fo the current time step
    if any(j==nt)
        PX_sol = PX;
        TT = tl(j);
    end

end

% Save Results
solution.T=TT;
solution.x=x;
solution.Xgrid=Xgrid;
solution.PTX=PX_sol;

end