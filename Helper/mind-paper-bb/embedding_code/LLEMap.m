classdef LLEMap < handle
    
    properties
        Xtrain; % (matrix)
            % Input training points (points x dimensions)
        Ytrain; % (matrix)
            % Output training points (points x dimensions)
        k; % (integer)
            % Number of nearest neighbors
        lambda; % (real)
            % Regularization parameter
        all_k % (vector)
            % All possible choices of k
        all_lambda % (vector)
            % All possible choices of lambda
        nfolds; % (integer)
            % Number of cross validation folds
    end % properties
    
    methods
        
        % Constructor
        function [ obj ] = LLEMap( k, lambda, nfolds )
            % Inputs:
            %   k
            %     Number of nearest neighbors. If a vector of values is
            %     given, one will be chosen by cross validation.
            %   lambda (Optional)
            %     Regularization parameter. If a vector of values is given,
            %     one will be chosen by cross validation. [Default: 0]
            %   nfolds (Optional)
            %     Number of cross validation folds. [Default: 10]
            
            obj.all_k = k(:)';
            
            if nargin < 2 || isempty(lambda)
                obj.all_lambda = 0;
            else
                obj.all_lambda = lambda(:)';
            end
            
            if nargin < 3 || isempty(nfolds)
                obj.nfolds = 10;
            else
                obj.nfolds = nfolds;
            end
        end % constructor
        
        
        % Fit mapping
        function [] = fit( obj, X, Y )
            % Inputs:
            %   X
            %     Matrix of training input points (points x dimensions)
            %   Y
            %     Matrix of training output points (points x dimensions)
            % Posteffects:
            %   Stores X and Y as object properties. Chooses k and lambda
            %   using cross validation if multiple values were provided.
            
            % store training data
            npts = size(X, 1);
            if size(Y, 1) ~= npts
                error([ ...
                    'Inputs and outputs must contain same number of ' ...
                    'points' ...
                ]);
            end
            obj.Xtrain = X;
            obj.Ytrain = Y;
            
            % set k and lambda if no cross validation requested
            nk = numel(obj.all_k);
            nlambda = numel(obj.all_lambda);
            if nk == 1 && nlambda == 1
                obj.k = obj.all_k;
                obj.lambda = obj.all_lambda;
                return;
            end
                
            % otherwise choose k and lambda by cross validation
            
            % construct cross validation folds
            cvp = cvpartition(npts, 'kfold', obj.nfolds);
            
            % find nearest neighbors for each cross validation fold
            nb = cell(1, obj.nfolds);
            D = squareform(pdist(X));
            kmax = max(obj.all_k);
            for it = 1 : obj.nfolds
                my_D = D(cvp.test(it), cvp.training(it));
                [~, idx] = sort(my_D, 2);
                nb{it} = idx(:, 1:kmax);
            end
            % nb{k}(i,j) = index of jth nearest neighbor in training
            % set of ith point in test set, for cross validation fold k
            
            % best values so far
            best_k = nan;
            best_lambda = nan;
            best_err = inf;

            % perform grid search
            ctr = 0;
            t = tic;
            
            % loop thru lambda
            for it_lambda = 1 : nlambda
                my_lambda = obj.all_lambda(it_lambda);

                % loop thru k
                for it_k = 1 : nk
                    my_k = obj.all_k(it_k);
                    
                    % squared test set prediction error
                    err = 0;
                    
                    % loop thru cross validation folds
                    for it_cv = 1 : obj.nfolds

                        % print status
                        ctr = ctr + 1;
                        if toc(t) >= 1
                            fprintf( ...
                                '%d/%d\n', ...
                                ctr, nlambda * nk * obj.nfolds ...
                            );
                            t = tic;
                        end
                        
                        % predict test set outputs
                        Yhat = obj.lle_mapping( ...
                            X(cvp.training(it_cv), :), ...
                            Y(cvp.training(it_cv), :), ...
                            X(cvp.test(it_cv), :), ...
                            nb{it_cv}(:, 1:my_k), ...
                            my_lambda ...
                        );
                        
                        % update error
                        Ytest = Y(cvp.test(it_cv), :);
                        err = err + sum((Ytest(:) - Yhat(:)).^2);
                    end % loop thru cross validation folds
                    
                    % keep parameters if they improve cost
                    if err < best_err
                        best_k = my_k;
                        best_lambda = my_lambda;
                        best_err = err;
                    end
                end % loop thru k
            end % loop thru lambda
            
            % store best parameters
            obj.k = best_k;
            obj.lambda = best_lambda;
        end % fit()
        
        
        % Map new input points into output space using LLE regression
        % method
        function [ Yi ] = lle_mapping( obj, X, Y, Xi, nb, lambda )
            % Inputs:
            %   X
            %     Training inputs
            %   Y
            %     Training outputs
            %   Xi
            %     New inputs
            %   nb
            %     Neighborhood matrix. nb(i,j) contains the index of the
            %     jth nearest training input to the ith new input
            %   lambda
            %     Regularization parameter (see lle_weights())
            % Outputs:
            %   Yi
            %     New points mapped into output space
            % Notes:
            %   All sets of points represented as matrics. Rows correspond
            %   to points, columns to dimensions.

            % data size
            npts_train = size(X, 1);
            npts_oos = size(Xi, 1);
            ndims_in = size(X, 2);
            ndims_out = size(Y, 2);

            % out-of-sample points mapped to output space
            Yi = zeros(npts_oos, ndims_out);

            % loop thru out-of-sample points
            for it = 1 : npts_oos

                % current out-of-sample point, neighboring training inputs
                my_Xi = Xi(it, :);
                my_nb = nb(it, :);

                % construct out-of-sample point in input space as a linear
                % combination of neighboring training points
                wts = obj.lle_weights(my_Xi, X(my_nb, :), lambda);

                % construct out-of-sample point in output space
                % apply same weights to image of neighbors in output space
                Yi(it, :) = wts(:)' * Y(my_nb, :);
            end
        end % lle_mapping()
    
    
        % Find weights for reconstructing a point as a linear combintion of
        % its neighbors, subject to constraint that weights sum to 1
        function [ wts ] = lle_weights( obj, x0, Xnb, lambda )
            % Inputs:
            %   x0
            %     The point to reconstruct. Given as a row vector.
            %   Xnb
            %     Neighboring points from which to reconstruct x0. Rows
            %     correspond to points, columns to dimensions.
            %   lambda
            %     Regularization term. Penalizes the summed squares of the
            %     weights (i.e. shrinks weights toward 0). Set to a small
            %     number (less than 1) if the problem is ill-conditioned
            %     (e.g. more points than dimensions).

            % number of neighbors
            n = size(Xnb, 1);

            % no need to continue if only a single neighbor given
            if n == 1
                wts = 1;
                return;
            end

            % calculate Gram matrix
            tmp = bsxfun(@minus, x0, Xnb);
            G = tmp * tmp';

            % condition Gram matrix to regularize weights
            if lambda > 0
                c = lambda^2 / n * trace(G);
                G = G + diag(ones(n,1) * c);
            end

            % solve for weight
            if det(G)<1e10
                try
                    wts = pinv(G) * ones(n, 1);          % with large entries in G, det tends to diverge and svd in turn does not converge.
                
                catch
                    disp('pinv failed, trying normal')
                    wts = inv(G) * ones(n, 1);
                end
            else
                try
                    scale0 = double(mean(G(:)));
                    fun = @(x) double(log(det(G./x)).^2);        % Find some scaling factor that works such that det(G/x) is not too bad...
                    scale = 3*fminbnd(fun, 1, 30*scale0);
                    wts = (pinv(G./scale)./scale) * ones(n, 1);
                catch
                    disp('Trying alternative...')
                    scale0 = nanmean(eig(G)); 
                    disp(['Found: nanmean(eig(G))=', num2str(scale0)])
                    fun = @(x) double(log(det(G./x)).^2);        % Find some scaling factor that works such that det(G/x) is not too bad...
                    disp(['Found: fun(1)=', num2str(fun(1))])
                    % make sure that the limits of the minimum function are
                    % finite values, and search in the broadest range
                    % possible
                    use_array = zeros(300,1);
                    testlims = linspace(1, 30*scale0, 300);
                    for t_idx = 1:length(testlims)
                        t_val = testlims(t_idx);
                        use_array(t_idx) = isfinite(fun(t_val));
                    end
                    [~,ind] = min(use_array); % Find the last value that is finite
                    if ind > 1
                        last_good_val = double(testlims(ind-1));
                    else
                        last_good_val = double(scale0);
                    end
                    disp(['Found: idx=', num2str(ind), ' val=', num2str(last_good_val)])
                    scale = 3*fminbnd(fun, 1, last_good_val);
                    disp(['Found: scale=', num2str(scale)])
                    disp(G)
                    try
                        wts = (pinv(G./scale)./scale) * ones(n, 1);
                    catch
                        wts = ones(n, 1);
                        disp("had to be cought")
                    end
                    disp(['Found: wts(1,1)=', num2str(wts(1,1)), ' ...continuing!'])
                end
            end
            if sum(abs(wts)) > 0
                wts = wts / sum(wts);
            else
                wts = ones(length(wts),1) / length(wts);
            end
        end % lle_weights()
    
    
        % Map new input points into output space 
        function [ Y ] = transform( obj, X )

            % check inputs
            if isempty(obj.Xtrain)
                error('Model must be fit first');
            end
            if size(X, 2) ~= size(obj.Xtrain, 2)
                error('Dimensionality of points doesn''t match model');
            end

            % for each point in X, find k nearest training inputs
            nb = knnsearch(obj.Xtrain, X, 'k', obj.k);

            % find mapped points
            Y = obj.lle_mapping(obj.Xtrain, obj.Ytrain, X, nb, obj.lambda);
        end % transform()
        
    end % methods
end % classdef
