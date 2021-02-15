function [dat, embed] = embedAsFunction(dat, parameters)
    
    % --------------------------------------------------------------------
    % parameters
    % --------------------------------------------------------------------
    assert(isfield(parameters.embed,'d'));
    assert(isfield(parameters.embed,'mode'));
    assert(strcmp(parameters.embed.mode, 'mds') || ...
        strcmp(parameters.embed.mode, 'cmds') || ...
        strcmp(parameters.embed.mode, 'hierarchical'));
    assert(~strcmp(dat.rwd.type, 'continuous'));
    if isfield(parameters.embed,'opts')
        opts = parameters.embed.opts;
    else
        opts = statset();
    end
    
    % --------------------------------------------------------------------
    % do not delete previous embedding info
    % --------------------------------------------------------------------
    if isfield(dat,'embed')
        if ~isfield(dat,'oldembed')
            dat.oldembed = {dat.embed};
        else
            dat.oldembed{end+1} = dat.embed;
        end
    end
    
    % --------------------------------------------------------------------
    % run mds on geodesic random walk distances
    % --------------------------------------------------------------------
    fprintf('\nsolving for mds embedding\n');
    dembed = parameters.embed.d;
    embed(numel(dembed)) = struct();
    switch parameters.embed.mode
        
        case 'mds' % non-classic multidimensional scaling
            for idim = 1:numel(embed)
                embed(idim).d = dembed(idim);
                embed(idim).opts = opts;
                embed(idim).mode = parameters.embed.mode;
            end
            yembed = cell(numel(embed),1);
            for idim = 1:numel(embed)
                fprintf('computing embedding %d of %d (d = %d)\n',...
                    idim,numel(embed),parameters.embed.d(idim));
                yembed{idim} = mdscale(dat.rwd.Dg, embed(idim).d,...
                    'criterion', 'sammon', 'Options', opts);
            end
            for idim = 1:numel(embed)
                embed(idim).y = yembed{idim};
            end
            
        case 'cmds' % classic multidimensional scaling
            for idim = 1:numel(embed)
                embed(idim).d = dembed(idim);
                embed(idim).opts = opts;
                embed(idim).mode = parameters.embed.mode;
            end
            yembed = cell(numel(embed),1);
            for idim = 1:numel(embed)
                fprintf('computing embedding %d of %d (d = %d)\n',...
                    idim,numel(embed),parameters.embed.d(idim))
                yembed{idim} = cmdscale(dat.rwd.Dg, embed(idim).d);
            end
            for idim = 1:numel(embed)
                embed(idim).y = yembed{idim};
            end
            
        case 'hierarchical' % hierarchical mode
            assert(isequal(dembed,min(dembed):max(dembed)));
            assert(isequal(dembed,sort(dembed)));
            for idim = numel(embed):-1:1
                fprintf('computing embedding %d of %d (d = %d)\n',...
                    idim,numel(embed),parameters.embed.d(idim))
                embed(idim).d = dembed(idim);
                embed(idim).opts = opts;
                embed(idim).mode = parameters.embed.mode;
                
                if idim == numel(embed)
                    % for largest dimension, use default initial conditons
                    embed(idim).y = mdscale(dat.rwd.Dg, embed(idim).d,...
                        'criterion', 'sammon', 'Options', opts);
                else
                    % for subsequent dimensions, initialize with previous output
                    embed(idim).y = mdscale(dat.rwd.Dg, embed(idim).d,...
                        'criterion', 'sammon', 'Options', opts, ...
                        'Start', embed(idim+1).y(:,1:(end-1)));
                end
            end
    end
    
    % --------------------------------------------------------------------
    % learn kernel mappings
    % --------------------------------------------------------------------
    if parameters.learnmapping
        % learn mapping from global pca space to manifold
        % using n.p. regression on landmark points
        for idim = 1:numel(embed)
            fprintf('computing mapping for dimension %d of %d (d = %d)\n',...
                idim,numel(embed),parameters.embed.d(idim))
            switch parameters.mapping.mode
                case 'lle'
                    % -------------------------------------------------------------------
                    % lle
                    % -------------------------------------------------------------------
                    k = parameters.mapping.k;
                    % choices of k for k nearest neighbors
                    lambda = parameters.mapping.lambda;
                    % choices of regularization parameter
                    nfolds_lle = parameters.mapping.nfolds_lle;
                    
                    % mapping from global pca space to embedding space
                    embed(idim).f2m.map = LLEMap(k, lambda, nfolds_lle);
                    embed(idim).f2m.map.fit(dat.lm.f, embed(idim).y);   %%% this is the bit.
                    
                    % map all points onto manifold
                    embed(idim).f2m.y = embed(idim).f2m.map.transform(dat.pca.f);
                    
                    % learn mapping from manifold to global pca space
                    % using lle regression method on landmark points
                    
                    % mapping from embedding space to global pca space
                    embed(idim).m2f.map = LLEMap(k, lambda, nfolds_lle);
                    embed(idim).m2f.map.fit(embed(idim).y, dat.lm.f);
                    
                    % reconstruct all points from estimated manifold coordinates
                    embed(idim).m2f.f = embed(idim).m2f.map.transform(embed(idim).f2m.y);
                    
                case 'wn'
                    % -------------------------------------------------------------------
                    % w-n kernel mode
                    % -------------------------------------------------------------------
                    h = parameters.mapping.h;
                    % choices of h
                    
                    % mapping from global pca space to embedding space
                    embed(idim).f2m.map = LLEMap_quick( h );
                    embed(idim).f2m.map.fit(dat.lm.f, embed(idim).y);
                    
                    % map all points onto manifold
                    embed(idim).f2m.y = embed(idim).f2m.map.transform(dat.pca.f);
                    
                    % learn mapping from manifold to global pca space
                    % using lle regression method on landmark points
                    
                    % mapping from embedding space to global pca space
                    embed(idim).m2f.map = LLEMap_quick( h );
                    embed(idim).m2f.map.fit(embed(idim).y, dat.lm.f);
                    
                    % reconstruct all points from estimated manifold coordinates
                    embed(idim).m2f.f = embed(idim).m2f.map.transform(embed(idim).f2m.y);
                    
                case 'tlle'
                    % -------------------------------------------------------------------
                    % lle with time-blocked cv
                    % -------------------------------------------------------------------
                    k = parameters.mapping.k;
                    % choices of k for k nearest neighbors
                    lambda = parameters.mapping.lambda;
                    % choices of regularization parameter
                    nfolds_lle = parameters.mapping.nfolds_lle;
                    
                    % mapping from global pca space to embedding space
                    t = dat.data.t(dat.lm.idx);
                    tmargin_sec = parameters.tmargin_sec;
                    embed(idim).f2m.map = LLEMap_time_block_cv(...
                        k, lambda, nfolds_lle);
                    embed(idim).f2m.map.fit(dat.lm.f, embed(idim).y, t, tmargin_sec);
                    
                    % map all points onto manifold
                    embed(idim).f2m.y = embed(idim).f2m.map.transform(dat.pca.f);
                    
                    % learn mapping from manifold to global pca space
                    % using lle regression method on landmark points
                    
                    % mapping from embedding space to global pca space
                    t = dat.data.t(dat.lm.idx);
                    tmargin_sec = parameters.tmargin_sec;
                    embed(idim).m2f.map = LLEMap_time_block_cv(...
                        k, lambda, nfolds_lle);
                    embed(idim).m2f.map.fit(embed(idim).y, dat.lm.f, t, tmargin_sec);
                    
                    % reconstruct all points from estimated manifold coordinates
                    embed(idim).m2f.f = embed(idim).m2f.map.transform(embed(idim).f2m.y);
                    
                case 'gp'
                    % -------------------------------------------------------------------
                    % gp
                    % -------------------------------------------------------------------
                    
                    if isfield(parameters.mapping, 'gp_options')
                        options = parameters.mapping.gp_options;
                    else
                        options = [];
                    end
                    
                    % mapping from global pca space to embedding space
                    embed(idim).f2m.map = LLEMap_gp();
                    embed(idim).f2m.map.fit(dat.lm.f, embed(idim).y, options);
                    
                    % map all points onto manifold
                    embed(idim).f2m.y = embed(idim).f2m.map.transform(dat.pca.f);
                    
                    % learn mapping from manifold to global pca space
                    % using lle regression method on landmark points
                    
                    % mapping from embedding space to global pca space
                    embed(idim).m2f.map = LLEMap_gp();
                    embed(idim).m2f.map.fit(embed(idim).y, dat.lm.f, options);
                    
                    % reconstruct all points from estimated manifold coordinates
                    embed(idim).m2f.f = embed(idim).m2f.map.transform(embed(idim).f2m.y);
            end
            % end switch
        end
        % end dimension loop
    end
    dat.embed = embed;
end % end embedAsFunction()
