%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot results of the PRESSTO algorithm:
% 1) Cost function evolution
% 2) Inducer profiles (stairs)
% 3) Final joint probability distribution
% 4) Final marginal distributions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PLOT(input_path)
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
    set(0, 'DefaultColorbarTickLabelInterpreter', 'latex');

    % ---- Folder mode ----
    if ~exist(input_path, 'dir')
        output_folder = fullfile(pwd, 'SIMULATIONS', input_path);
    else
        output_folder = input_path;
    end
    
    results_file = fullfile(output_folder, 'Results.mat');
    config_file  = fullfile(output_folder, 'config.mat');

    if ~isfile(results_file)
        error('PLOT:MissingResults', ...
            'Results.mat not found in folder: %s', output_folder);
    end
    if ~isfile(config_file)
        error('PLOT:MissingConfig', ...
            'config.mat not found in folder: %s', output_folder);
    end

    % Load results and configuration
    load(results_file, 'Results');
    load(config_file, 'config');
    meshxt = config.meshxt;

    assert(isfield(Results,'PX'), ...
        'Results.PX not found. Invalid Results struct.');

    %% ============================= PLOTTING =================================
    %% ------------------ TIME ------------------ %%
    t = Results.t;

    %% ------------------ COST FUNCTION ------------------ %%
    figure;
    plot(t, Results.cost, 'LineWidth', 2);
    grid on;
    xlabel('Time');
    ylabel('Cost function');
    title('Cost function evolution');
    set(gca, 'FontSize', 16);

    %% ------------------ INDUCER PROFILES ------------------ %%
    u = Results.u;
    n_gene = size(u,2);
    u_plot = [u; u(end, :)]; 
    
    figure('Color', 'w');
    for i = 1:n_gene
        subplot(n_gene, 1, i)
        stairs(t, u_plot(:, i), 'LineWidth', 2);
        grid on;
        set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');
        ylabel(['Inducer ' num2str(i)], 'Interpreter', 'latex');
        xlim([t(1), t(end)]); 
        
        if i == 1
            title('Inducer profiles', 'Interpreter', 'latex', 'FontSize', 16);
        end
        if i == n_gene
            xlabel('Time', 'Interpreter', 'latex');
        end
    end

    %% ------------------ DISTRIBUTIONS ------------------ %%
    PX = Results.PX{end};
    dim = sum(size(PX) > 1);
    n_gene = size(meshxt.Prot_mesh, 1);

    % Build spatial grid
    if gpuDeviceCount > 0 
        meshxt.Prot_mesh_gpu = gpuArray(meshxt.Prot_mesh);
        iN = cell(n_gene,1); 
        x = cell(n_gene,1);
        Xgrid = cell(n_gene,1);
        for i = 1:n_gene
            iN{i} = meshxt.Prot_mesh_gpu(i,3) + 1;
            x{i} = linspace(meshxt.Prot_mesh_gpu(i,1), meshxt.Prot_mesh_gpu(i,2), iN{i}); 
        end
        [Xgrid{1:n_gene}] = ndgrid(x{:});
    else 
        iN = cell(n_gene,1);
        x = cell(n_gene,1);
        for i = 1:n_gene
            iN{i} = meshxt.Prot_mesh(i,3) + 1;
            x{i} = linspace(meshxt.Prot_mesh(i,1), meshxt.Prot_mesh(i,2), iN{i});
        end
        Xgrid = cell(n_gene,1);
        [Xgrid{1:n_gene}] = ndgrid(x{1:n_gene});
    end

    %% ---- Joint distribution ---- %%
    if dim == 1
        figure;
        plot(x{1}, PX, 'LineWidth', 2);
        set(gca, 'FontSize', 16);
        xlabel('Protein');
        ylabel('Probability');
        title(sprintf('Final Controlled Distribution at $t = %g$', meshxt.Time_mesh(2)));
        grid on;
    elseif dim == 2
        figure;
        mesh(Xgrid{1}, Xgrid{2}, PX); 
        set(gca, 'FontSize', 16);
        xlabel('Protein 1');
        ylabel('Protein 2');
        zlabel('Probability');
        title(sprintf('Final Controlled Distribution at $t = %g$', meshxt.Time_mesh(2)));
        shading interp;
    elseif dim == 3
        % --- Select isosurface levels ---
        PXmax = max(PX(:));
        alpha_min = 0.1;
        alpha_max = 0.99;
        nlevels = 9;
        levels = linspace(alpha_min*PXmax, alpha_max*PXmax, nlevels);
        
        % --- Colormap ---
        cmap = parula(numel(levels));
        
        figure; hold on;
        
        % FaceAlpha proporcional a los niveles
        face_alpha = linspace(alpha_min, alpha_max, nlevels);
        
        for i = 1:nlevels
            p = patch(isosurface(Xgrid{1}, Xgrid{2}, Xgrid{3}, PX, levels(i)));
            p.FaceColor = cmap(i,:);
            p.EdgeColor = 'none';
            p.FaceAlpha = face_alpha(i);
        end
        % --- Visualization settings ---
        camlight;
        lighting gouraud;
        axis tight;
        axis equal;
        grid on;
        set(gca, 'FontSize', 16);
        xlabel('Protein 1');
        ylabel('Protein 2');
        zlabel('Protein 3');
        title(sprintf('Final Controlled Distribution at $t = %g$', meshxt.Time_mesh(2)));
    end

    %% ---- Marginal distributions ---- %%
    if dim > 1
        for i = 1:dim
            PX_i = PX;
            for j = 1:dim
                if j ~= i
                    PX_i = trapz(x{j}, PX_i, j);
                end
            end
            PX_i = squeeze(PX_i);
            PX_i = PX_i / trapz(x{i}, PX_i);
            figure;
            plot(x{i}, PX_i, 'LineWidth', 2);
            set(gca, 'FontSize', 16);
            xlabel(['Protein ' num2str(i)]);
            ylabel('Probability');
            title(['Final marginal distribution of protein ' num2str(i)]);
            grid on;
        end
    end
end