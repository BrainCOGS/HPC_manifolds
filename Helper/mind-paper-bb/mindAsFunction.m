function dat = mindAsFunction(data,parameters)
    
    % saves output in struct 'dat' with fields/sub-structs:
    %   data: original data
    %   pca: global pca
    %   forest: ppca forest for estimating transition probabilities
    %   lm: landmark points
    %   tp: transition probabilities
    %   rwd: random walk distances
    
    if ~isfield(parameters,'normalizepTF')
        parameters.normalizepTF = true;
    end
    
    % all analysis outputs will be stored in this struct
    dat = struct;
    
    % defined in environment:
    %   t: time of each sample
    %   x: trajectory in environment (time x dimensions)
    %   f: network activity over time (time x neurons)
    
    assert( size( data.t , 2 ) == 1);
    nt = size( data.t , 1 );
    assert( size( data.f , 1 ) == nt );
    
    dat.data.t = data.t;
    dat.data.f = data.f;
    
    % select contiguous samples, according to selected time step
    dat.dt = parameters.dt;
    dt = dat.dt;
    assert(dt == 1 || dt == 0); % do want to eventually make it consistent across dt > 1
    if dt > 0
        % select contiguous samples
        delta_t = diff(dat.data.t);
        si = median(delta_t);
        dat.data.tsel = find(delta_t < 1.5*si);
        % this will work for our current setup. Eventually just want to
        % include a user-input mechanism for determining "1.5*si"
    else
        dat.data.tsel = 1:numel(dat.data.t);
    end
    
    % preprocess data with global pca
    
    % --- parameters ---
    % dat.pca.n = size(f,2);
    dat.pca.n = parameters.pca.n;
    % how many principal components to retain (or fraction of variance)
    % ------------------
    
    % run pca
    dat.pca.model = MyPCA;
    dat.pca.model.fit(dat.data.f);
    if ~isnan( dat.pca.n )
        dat.pca.f = dat.pca.model.transform(dat.data.f, dat.pca.n);
    else
        dat.pca.f = dat.data.f;
    end    
    
    % select landmark points
    
    % --- parameters ---
    dat.lm.n = parameters.lm.n;
    % how many landmark points
    % ------------------

    %%% Edit by Manuel: Replaced landmarks code with the greedy algoritm
%         % initialize with greedy vector quantization
%         as = greedyvq(dat.pca.f, dat.lm.n);
%         % prune duplicate landmarks
%         % fine tune with k medoids
%         [~, dat.lm.idx] = kmedoids_d( ...
%             squareform(pdist(dat.pca.f)), ...
%             dat.lm.n, ...
%             'init', as ...
%             );
    
    dat.lm.idx = get_landmarks(dat.pca.f,dat.lm.n);
    dat.lm.idx = sort(dat.lm.idx); % sort landmarks so they are in time-increasing order
    % Note that for lmf = 1, the last point cannot be a landmark, therefore
    % point "1" is duplicated.
   
    % prune landmarks so they are separated in time, if requested
    if ~isempty(parameters.prune_lm_by_time) && ~isnan(parameters.prune_lm_by_time)
        tlm = dat.data.t( dat.lm.idx );
        tlm_good = tlm;
        while min(diff(tlm_good)) < parameters.prune_lm_by_time
            [~,is] = min(diff(tlm_good));
            is = is + round(rand());
            tlm_good(is) = [];
        end
        dat.lm.idx = dat.lm.idx(ismember(tlm, tlm_good));
    end
    fprintf('Starting at %d in the parameters, \n', parameters.lm.n);
    fprintf('There are %d landmark points after all pruning\n',numel(dat.lm.idx));
    dat.lm.n = numel(dat.lm.idx);

    parameters.lm.n = dat.lm.n;
    % idx contains indices of points to be used as landmarks
    
    % pca projections and environment location for each landmark point
    dat.lm.f = dat.pca.f(dat.lm.idx, :);
        
    % --------------------------------------------------------------------
    % train forest
    % --------------------------------------------------------------------    
    % initialize ppca forest
    node = HPPCARegressionNode( ...
        'dim_criterion', parameters.dim_criterion, ...
        'ndir', parameters.ndir, ...
        'min_leaf_pts', parameters.min_leaf_pts ...
        );
    
    dat.forest = DForest( ...
        node, ...
        'ntrees', parameters.ntrees , ...
        'verbose', parameters.verbose ...
        );
    % train forest to predict next state as a function current state, next = x0 + dt    
    trees = dat.forest.trees;
    ntrees = numel(trees);
    it = 1;
    while it <= ntrees
        fprintf('training tree %d/%d\n', it, ntrees);
        try %sometimes svd does not converge. if this happens, just try again.
        trees(it).train( ...
            dat.pca.f(dat.data.tsel, :), ...
            dat.pca.f(dat.data.tsel + dt, :) ...
            );
        trees(it).patience = dat.forest.patience;
        it = it + 1;
        catch
            fprintf('tree %d/%d. FAILED. Redoing...\n', it, ntrees);
            trees(it) = DTree(node, 'patience', dat.forest.patience);  %reinitialize the tree.
        end
    end    
    dat.forest.ndims_in = size(dat.pca.f(dat.data.tsel, :), 2);
    
    % --------------------------------------------------------------------
    % compute transition probabilities between landmark points
    % --------------------------------------------------------------------
    % log transition probabilities for each tree
    
    if dat.lm.n < 3000
        % This is faster, but doesn't work with many landmarks, because
        % matrix becomes too big.
        lp_all = zeros(dat.lm.n, dat.lm.n, dat.forest.ntrees);
        for it = 1 : ntrees
            fprintf('computing log transition probabilities for tree %d/%d\n', ...
                it, ntrees);
            lp_all(:, :, it) = trees(it).lpygx(dat.lm.f, dat.lm.f, true);
        end
        dat.forest.trees = trees;
        dat.tp.lp_all = lp_all;    
        % lp_all(i,j,k) = log probability density that the next network state is
        % jth landmark point, given that the current state is ith landmark point,
        % according to tree k
        log_p_matrix = median(dat.tp.lp_all,3);
    else
        log_p_matrix = zeros(dat.lm.n,dat.lm.n);
        for x_idx = 1:dat.lm.n
            fprintf('computing log transition probabilities %d/%d\n',x_idx, dat.lm.n);
            xx = dat.lm.f(x_idx,:);
            lp_all = zeros(ntrees,dat.lm.n);
            for it = 1 : ntrees
                lp_all(it,:) = trees(it).lpygx(xx, dat.lm.f, true);
            end
            log_p_matrix(x_idx,:) = median(lp_all,1);
        end
        dat.forest.trees = trees;
    end
    
    %%% Code from manuel inserted for normalization and log-p treatment
    dat.tp.lp = zeros(size(log_p_matrix));
    for i = 1:dat.lm.n
        dat.tp.lp(i,:) =  log_p_matrix(i,:) - log(sum(exp(log_p_matrix(i,:))));
    end
    
    ma = max(dat.tp.lp(:));
    f1 = sum(dat.tp.lp(:) >= (ma-10) );
    f2 = sum((dat.tp.lp(:) >= (ma-20)).*dat.tp.lp(:) < (ma-10) );
    disp("F1/F2 ratio - raw: " + num2str(f1/f2))
    
    dat.tp.f1f2_raw = f1/f2;
    
    if f1/f2 > 1.5  % "1.5" is heuristic. For good embeddings, f1/f2 <~1, 
        disp("Likelihoods unstable. Trying a clip.")

        % Make a conditioned matix without the infs:
        log_p_conditioned = log_p_matrix;
        max_allowed = log(realmax)/2;      % This is a matlab constrained. Make sure sum(exp(...)) remains finite for a big number of neurons
        log_p_conditioned(log_p_conditioned>max_allowed) = max_allowed;
        dat.tp.lp = zeros(size(log_p_matrix));
        for i = 1:dat.lm.n
            dat.tp.lp(i,:) =  log_p_conditioned(i,:) - log(sum(exp(log_p_conditioned(i,:))));
        end
        % Test, after this normalization "sum(exp(dat.tp.lp),2)" has to be one.
        
        ma = max(dat.tp.lp(:));
        f1 = sum(dat.tp.lp(:) >= (ma-10) );
        f2 = sum((dat.tp.lp(:) >= (ma-20)).*dat.tp.lp(:) < (ma-10) );
        disp("F1/F2 ratio - clipped: " + num2str(f1/f2));
        dat.tp.f1f2_clipped = f1/f2;
    end
    
    if f1/f2 > 1.5
        disp("Likelihoods still unstable. Trying a clip + shift.")
        dat.tp.lp = zeros(size(log_p_matrix));
        log_p_conditioned = log_p_matrix;
        log_p_conditioned = log_p_conditioned - (max(log_p_matrix(:)) - max_allowed);
        lptest = zeros(size(log_p_matrix));
        for i = 1:dat.lm.n
            dat.tp.lp(i,:) =  log_p_conditioned(i,:) - log(sum(exp(log_p_conditioned(i,:))));
        end
        ma = max(dat.tp.lp(:));
        f1 = sum(dat.tp.lp(:) >= (ma-10) );
        f2 = sum((dat.tp.lp(:) >= (ma-20)).*dat.tp.lp(:) < (ma-10) );
        disp("F1/F2 ratio - clipped+shifted: " + num2str(f1/f2));
        
        dat.tp.f1f2_clipped_shifted = f1/f2;
    end
    
    
%%%% This is S+R original:
%     maxval = max(dat.tp.lp_all(:));
%     dat.tp.lp = log(median(exp(dat.tp.lp_all - maxval), 3)) + maxval;
%     % lp(i,j) = log of median i->j transition probability over trees
%     % compute this way to avoid floating point underflow
%     % normalized transition probabilities
%     maxval = max(dat.tp.lp(:));
%     lpsum = log(sum(exp(dat.tp.lp - maxval), 2)) + maxval;
%     % lpsum(i) = log sum_{j=1:nc} i->j transition probability
%     dat.tp.np = exp(bsxfun(@minus, dat.tp.lp, lpsum));
%     % np(i,j) = log of normalized i->j transition probability
%     % normalization is such that transition probabilities sum to 1 over
%     % discrete set of states being considered

%%%% This is one of S+R testroutines.
%
%     % remove disconnected landmark points
%     % there should only be a few scattered ones. more indicates a disconnected
%     % manifold or a problem (e.g. landmark points to far away from each other)
%     
%     % find connected components of transition probability graph
%     [~, as] = graphconncomp(sparse(dat.tp.np));
%     
%     % how many landmark points in each connected component
%     count = accumarray(as(:), 1);
%     
%     % there are disconnected landmark points
%     % if graph contains multiple connected components
%     % if so, remove landmark points
%     if numel(count) > 1
%         
%         % keep landmark points in the largest connected component
%         [~, ccmax] = max(count);
%         sel = as == ccmax;
%         
%         % remove the others
%         nbad = nnz(~sel);
%         dat.lm.n = dat.lm.n - nbad;
%         dat.lm.idx = dat.lm.idx(sel);
%         dat.lm.f = dat.lm.f(sel, :);
%         dat.tp.lp_all = dat.tp.lp_all(sel, sel, :);
%         dat.tp.lp = dat.tp.lp(sel, sel);
%         dat.tp.np = dat.tp.np(sel, sel);
%         
%         % renormalize transition probabilities
%         dat.tp.np = bsxfun(@rdivide, dat.tp.np, sum(dat.tp.np, 2));
%         
%         fprintf( ...
%             'Removed %d landmark points (%d connected components)\n', ...
%             nbad, numel(count) - 1 ...
%             );
%     else
%         fprintf('No landmark points removed\n');
%     end        
    

    % --- parameters ---
    dat.rwd.type = parameters.rwd.type;
    % type of random walk ('continuous' or 'discrete')
    dat.rwd.sym = parameters.rwd.sym;
    % how to symmetrize global distances ('avg' or 'min')
    dat.rwd.all_geo = parameters.rwd.all_geo;
    % if true, all distances will be geodesic distances. otherwise,
    % distances between connected points will be ideal local distances and
    % distances between non-connected points will be filled in with
    % geodesic distances.
    % ------------------
    
    switch dat.rwd.type
        
        % continuous random walk
        case 'continuous'
            dat.rwd.d = parameters.rwd.d;
            % dimensionality of space in which random walk is performed
            % (only used for continuous random walk)
            
            dat.rwd.var_scale = parameters.rwd.var_scale;
            % variance of diffusion kernel, expressed as fraction of maximum
            % possible variance (only used for continuous random walk; shouldn't
            % matter much)
            
            % log transition probability densities
            %L = log(dat.tp.np);
            L = dat.tp.lp;
            
            % variance of diffusion kernel
            Pmax = max(exp(L(:)));
            sigma2_max = 1/(2*pi*Pmax^(2/dat.rwd.d));
            sigma2 = sigma2_max * dat.rwd.var_scale;
            
            % ideal local distances
            dat.rwd.Dl = sqrt( ...
                -2 * sigma2 * L ...
                - dat.rwd.d * sigma2 * log(2*pi*sigma2) ...
                );
            dat.rwd.Dl(1 : size(L, 1)+1 : end) = 0;
            
            % discrete random walk
        case 'discrete'
            % transition probabilities (rows normalized to sum to 1)
            % ideal local distances
            
%  for Sam analysis: true! thresholds probabilities
%             if parameters.normalizepTF
%                 P = dat.tp.np;
%                 dat.rwd.Dl = sqrt(-log(P));
%             else
                P = dat.tp.lp;
                dat.rwd.Dl = sqrt(-P);
%             end

            dat.rwd.Dl(1 : size(P, 1)+1 : end) = 0;
            
        otherwise
            error('Random walk type must be ''continuous'' or ''discrete''');
    end
    
    % normalize by maximum distance for stability
    maxval = max(dat.rwd.Dl(isfinite(dat.rwd.Dl)));
    dat.rwd.Dl = dat.rwd.Dl ./ maxval;
    
    % geodesic distances (global)
    dat.rwd.Dl = real(dat.rwd.Dl); % why do these things get to be complex? ...
    A = isfinite(dat.rwd.Dl);
    % adjacency matrix; A(i,j) = 1 if landmark points i,j are connected    
    dat.rwd.Dg = graphallshortestpaths( ...
        sparse(A), ...
        'weights', dat.rwd.Dl(A) ...
        ); % geodesic distance matrix
    fprintf('There are %d lm points\n', ...
        size(dat.rwd.Dg,1));
    
    % symmetrize distances
    switch dat.rwd.sym
        case 'avg'
            dat.rwd.Dg = (dat.rwd.Dg + dat.rwd.Dg') / 2;
        case 'min'
            dat.rwd.Dg = min(dat.rwd.Dg, dat.rwd.Dg');
        otherwise
            error('Symmetrization method must be ''avg'' or ''min''');
    end
    
    %%% Edit by Manuel: if there are zero distances, these are removed from Dg.
    % This is necessary for datasets with very small probabilities that end
    % up after normalization as 1/0 and mess up the subsequent MDS.
    % this removes the landmarks and the entries from distance matrix.
    Dnew = dat.rwd.Dg;
    del_idx = [];
    for ind = find(dat.rwd.Dg==0)'
        if mod(ind,dat.lm.n) ~= ceil(ind/dat.lm.n)
            del_idx = [del_idx; mod(ind,dat.lm.n)];
            disp('Zero distance between two landmarks found. lm removed:')
            disp(mod(ind,dat.lm.n))
        end
    end
    bad_lms = del_idx(2:2:end); % always occurs in pairs, only one landmark has to go.
    bad_lms = [bad_lms; find(sum(isinf(Dnew), 2) >= (dat.lm.n-1) )]; %and the infinite distances.

    Dnew(bad_lms,:) = [];  
    Dnew(:,bad_lms) = [];
    dat.rwd.Dg = Dnew;
    % update landmarks
    sel = ~ismember(1:dat.lm.n, bad_lms);
    dat.lm.n = dat.lm.n - length(bad_lms);
    dat.lm.idx = dat.lm.idx(sel);
    dat.lm.f = dat.lm.f(sel, :);
    
    
    
    % if requested, use ideal local distances between conencted points
    if ~dat.rwd.all_geo
        if strcmp(dat.rwd.sym, 'avg')
            tmp = (dat.rwd.Dl + dat.rwd.Dl') / 2;
        elseif strcmp(sym_type, 'min')
            tmp = min(dat.rwd.Dl, dat.rwd.Dl');
        end
        dat.rwd.Dg(isfinite(tmp)) = tmp(isfinite(tmp));
    end    
    
end