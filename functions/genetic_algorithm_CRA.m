function [params] = genetic_algorithm_CRA(params, ga_params)

    mutation = ga_params.mutation;
    crossover = ga_params.crossover;
    
    pop_struct.val = nan(ga_params.n_genes, 1);
    pop_struct.cost = nan;
    old_pop_struct = repmat(pop_struct, ga_params.n_pop, 1);
    
    best_sol_struct.It = nan;
    best_sol_struct.cost = nan;
    best_sol_struct.time = nan;
    best_sol_struct.val = nan(ga_params.n_genes, 1);
    best_sol_struct = repmat(best_sol_struct, ga_params.iter_max+1, 1);
    
    % Initialize population
    tStart = tic;
    for i = 1:ga_params.n_pop
        chrom_uniform = linspace(0, 1, ga_params.n_genes);
    
        % add random gaussian noise
        sigma = ((chrom_uniform(2) - chrom_uniform(1))/5);
        chrom_curr = chrom_uniform + normrnd(0, sigma, [1, ga_params.n_genes]);
    
        chrom_curr = fix_boundaries(chrom_curr);
    
        old_pop_struct(i).val = chrom_curr;
    end
    
    new_pop_struct = repmat(pop_struct, ga_params.n_pop, 1);
    
    for i = 1:ga_params.iter_max
        ga_params.iter = i;
        tStart = tic;
        
        % evaluate pop cost
        parfor ii = 1:ga_params.n_pop
            old_pop_struct(ii).cost = cost(old_pop_struct(ii).val, params);
        end
        
        % penalty for non solved chromosomes
        best_sol_struct(i).It = i;
        [old_pop_struct(isnan([old_pop_struct.cost].')).cost] = deal(max([old_pop_struct.cost].')); 
        
        
        % offset + rank scaling
        old_pop_fit = -[old_pop_struct.cost].';
        old_pop_fit = old_pop_fit - min(old_pop_fit);
        old_pop_fit = rank_scaling(old_pop_fit);
    
        % elitism
        [~, idx] = sort([old_pop_struct.cost]);
        new_pop_struct(1:ga_params.n_elite) = old_pop_struct(idx(1:ga_params.n_elite));
    
        % selection + shuffle
        n_selections = 2*ga_params.n_crossover;
        sel_idx = stochastic_universal_selection(old_pop_fit, n_selections);
        sel_idx = sel_idx(randperm(length(sel_idx)));
        sel_idx = reshape(sel_idx, [], 2);
        
        % crossover
        crossover_out = [];
        for ii = 1:size(sel_idx, 1)
            idx1 = sel_idx(ii, 1);
            idx2 = sel_idx(ii, 2);
            
            chrom_0 = crossover(old_pop_struct(idx1).val, old_pop_struct(idx2).val);
            chrom_0 = fix_boundaries(chrom_0);
    
            tmp = pop_struct;
            tmp.val = chrom_0;
            crossover_out = [crossover_out; tmp];
        end
        new_pop_struct((ga_params.n_elite + 1):ga_params.n_pop) = crossover_out(1:(ga_params.n_pop - ga_params.n_elite));
            
        % mutation
        for ii = (ga_params.n_elite + 1):ga_params.n_pop
            if rand <= ga_params.p_mu
                chrom_0 = mutation(new_pop_struct(ii).val, ga_params);
                chrom_0 = fix_boundaries(chrom_0);
                
                new_pop_struct(ii).val = chrom_0;
            end
        end
    
    
        % save and disp results
        [~, idx] = min([new_pop_struct.cost]);
        best_sol_struct(i).val = chrom2Rm(new_pop_struct(idx).val, params)*1e2;
        best_sol_struct(i).cost = new_pop_struct(idx).cost;
        tEnd = toc(tStart);
        best_sol_struct(i).time = tEnd;
        T = struct2table(best_sol_struct, 'AsArray', true);
        clc; disp(T(1:i, :));
    
        % update next generation
        old_pop_struct = new_pop_struct;
        new_pop_struct = repmat(pop_struct, ga_params.n_pop, 1);
    end


    params.Rm = (best_sol_struct(ga_params.iter_max).val.')*1e-2;
    params.Nm = ceil(4*pi*params.Rm*max(params.fband)/params.c);
    params.Nm(params.Nm == 0) = 1;
    params = update_params(params);

    
end


function costval = cost(chrom, params)
    params.Rm = chrom2Rm(chrom, params).';
    params.Nm = ceil(4*pi*params.Rm*max(params.fband)/params.c);
    params.Nm(params.Nm == 0) = 1;
    params = update_params(params);
    [costval, ~] = calc_proposed_FIR_beamformer(params);
end


function chrom = Rm2chrom(Rm, params)
    chrom = Rm/params.Rmax;
end


function Rm = chrom2Rm(chrom, params)
    Rm = params.Rmin + chrom*(params.Rmax - params.Rmin);
end


function chrom = fix_boundaries(chrom)
    % clipping
    chrom(chrom >= 1) = 1;
    chrom(chrom <= 0) = 0;
    chrom = sort(chrom);
end