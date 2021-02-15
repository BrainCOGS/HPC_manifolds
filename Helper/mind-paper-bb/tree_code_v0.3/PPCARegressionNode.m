% PPCARegressionNode class
% Implements nodes of a regression tree with multivariate outputs. Uses
% hyperplanes to partition the input space. Fits a probabilistic PCA model
% to the outputs for each node.

% References:
%   Tipping and Bishop (1999). Probabilistic principal component analysis.
%   Minka (2000). Automatic choice of dimensionality for PCA.


classdef PPCARegressionNode < BSHNode
    
    properties
        
        % Local model parameters
        d; % (scalar)
            % Output dimensionality
        q; % (scalar)
            % PPCA dimensionality
        mu; % (1 x q vector)
            % Mean
        U; % (d x q matrix)
            % Top q eigenvectors of PPCA covariance matrix
        lambda; % (q x 1 vector)
            % Top q eigenvalues of PPCA covariance matrix
        sigma2; % (scalar)
            % PPCA noise variance
        ldt; % (scalar)
            % First two terms for the log density of a single data point
            % -d/2*log(2*pi) - 1/2*log(det(C))
            % where C is the PPCA covariance matrix
        
        % Diagnostic outputs
        eigs; % (d x 1 vector)
            % Eigenvalues of the sample covariance matrix
        evidence; % (vector)
            % If using Bayesian model selection to choose the
            % dimensionality, evidence(j) contains the log evidence for
            % dimensionality j. Otherwise empty.
            
        % Shared parameters stored in params object:
        %
        % dim_criterion
        %   How to choose dimensionality of local PPCA models. If >= 1,
        %   specifies a fixed dimensionality. If between 0 and 1,
        %   dimensionality will bechosen to account for at least this
        %   fraction of the variance. If <= 0, dimensionality will be
        %   chosen to maximize the Bayesian evidence.
        % split_evidence
        %   Determines loss function to use when optimizing splits. If
        %   true, uses average negative log Bayesian evidence. Otherwise
        %   uses average negative log likelihood.
        % ngrid
        %   Number of points for initial grid search when optimizing
        %   splitting threshold.
        % ndir
        %   Number of split directions to try when optimizing splits.
        % axis_parallel
        %   If true, splitting hyperplanes will be axis parallel. If false,
        %   they'll be randomly oriented.
        % random_split
        %   If true, a random threshold and direction will be chosen for
        %   each split. Otherwise, they'll be optimized.
        % max_depth
        %   Maximum depth of the tree (root node has depth 0)
        % min_leaf_pts;
        %   Minimum number of training points per leaf node
        % min_leaf_loss;
        %   A node can only be split if the loss of its local model exceeds
        %   this value.

        
    end % properties
    
    methods
        
        % Constructor
        function [ obj ] = PPCARegressionNode( varargin )
            % Inputs:
            %   If the empty matrix is given as input, creates an empty
            %   PPCARegressionNode.
            %
            %   Otherwise (if creating a new root node), parameters can be
            %   specified as optional name/value pairs:
            %     'dim_criterion'
            %       How to choose dimensionality of local PPCA models.
            %       If >= 1, specifies a fixed dimensionality. If between
            %       0 and 1, dimensionality will bechosen to account for
            %       at least this fraction of the variance. If <= 0,
            %       dimensionality will be chosen to maximize the Bayesian
            %       evidence. [Default: 0.95]
            %     'split_evidence'
            %       Determines loss function to use when optimizing splits.
            %       If true, uses average negative log Bayesian evidence
            %       (requires dim_criterion < 0). Otherwise uses average
            %       negative log likelihood. [Default: false]
            %     'ngrid'
            %       Number of points for initial grid search when
            %       optimizing splitting threshold.
            %     'ndir'
            %       Number of split directions to try when optimizing
            %       splits. [Default: 1]
            %     'axis_parallel'
            %       Perform axis-parallel splits if true, otherwise split
            %       in arbitrary directions. [Default: false]
            %     'random_split'
            %       If true, a random threshold and direction will be
            %       chosen for each split. Otherwise, they'll be optimized.
            %       [Default: false]
            %     'max_depth'
            %       Maximum depth of the tree (root node has depth 0)
            %       [Default: inf]
            %     'min_leaf_pts'
            %       Minimum number of training points per leaf node.
            %       [Default: 1]
            %     'min_leaf_loss'
            %       A node can only be split if the loss of its local model
            %       exceeds this value. [Default: -inf]
            %       Number of split directions to try. [Default: 1]
            % Outputs:
            %   Returns a PPCARegressionNode object whose shared parameter
            %   object ('params' property) has been set, if requested. All
            %   other properties are empty.
            
            % set shared parameters if not creating an empty node
            if ~(nargin == 1 && isempty(varargin{1}))
                default_args = struct( ...
                    'dim_criterion', 0.95, ...
                    'split_evidence', false, ...
                    'ngrid', Inf, ...
                    'ndir', 1, ...
                    'axis_parallel', false, ...                    
                    'random_split', false, ...
                    'max_depth', inf, ...
                    'min_leaf_pts', 1, ...
                    'min_leaf_loss', -inf ...
                );
                obj.params = HandleVar();
                obj.params.data = ...
                    parse_args(varargin, default_args, false);
            end
        end % constructor
        
        
        % Create a new, empty node with same shared parameters
        function [ obj2 ] = spawn( obj )
            % Outputs:
            %   obj2
            %     An empty node. Its 'params' HandleVar object is copied
            %     from obj.
            
            obj2 = PPCARegressionNode([]);
            obj2.params = obj.params;
        end
        
        
        % Train local PPCA model for the outputs
        function [ loss ] = train( obj, x, y )
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions.
            %   y
            %     Output data matrix. Rows correspond to points, columns to
            %     dimensions.
            %   mu (Optional)
            %     Mean of y. Will be computed if not given.
            %   eigval (Optional)
            %     Eigenvalues of the sample covariance matrix of y, given
            %     as a vector in descending order. Will be computed if
            %     not given.
            %   eigvec (Optional)
            %     Corresponding eigenvectors (stored on the columns).
            %     Will be computed if not given.
            % Outputs:
            %     Value of the loss function for the local model on the
            %     training data.
            
            % number of points, output dimensions
            [npts, nd] = size(y);
            if nd <= 1
                error('Outputs must be multivariate');
            end
            obj.d = nd;
            obj.train_pts = npts;
            
            % maximum likelihood PPCA mean is the sample mean
            % (T&B just below equation 5)
            mu = ones(1, npts) * y ./ npts;
            obj.mu = mu;
            
            % eigendecomposition of the sample covariance matrix
            % use eigendecomposition of the covariance matrix or svd of
            % the data matrix (whichever is faster)                        
            if npts > nd                               
                [eigvec, eigval] = eig(cov(y, 1), 'vector');
                [eigval, ord] = sort(eigval, 'descend');
                eigvec = eigvec(:, ord);
            else
                [eigvec, S] = svd(bsxfun(@minus, y, mu)', 0);
                eigval = diag(S) .^ 2 ./ npts;
                
            end
                % eigval(j) contains the jth largest eigenvalue
                % eigvec contains corresponding eigenvectors on the columns
            obj.eigs = eigval;
            
            % fit ppca model
            [q, obj.sigma2, obj.ldt, obj.evidence, loss] = ...
                obj.fit_ppca(npts, obj.mu, eigval, eigvec);
            obj.q = q;
            obj.train_loss = loss;
            
            % top q eigenvalues/eigenvectors of the PPCA covariance matrix
            obj.U = eigvec(:, 1:q);
            obj.lambda = eigval(1:q);
        end % train()
        
        
        % Fit PPCA model to training data
        function [ q, sigma2, ldt, evidence, loss ] = ...
                fit_ppca( obj, n, mu, eigval, eigvec )
            % Inputs:
            %   n
            %     Number of data points.
            %   mu
            %     Mean of data.
            %   eigval
            %     Eigenvalues of the sample covariance matrix, stored in
            %     a vector in descending order.
            %   eigvec
            %     Corresponding eigenvectors (stored on the columns).
            % Outputs:
            %   q
            %     Dimensionality of the latent variables. Chosen using
            %     method specified in obj.params
            %   sigma2
            %     Maximum likelihood noise variance.
            %   ldt
            %     Sum of the first two terms for the log density of a data
            %     point: -d/2*log(2*pi) - log(det(C))/2
            %   evidence
            %     If using Bayesian model selection to choose
            %     dimensionality, evidence(j) contains the log evidence for
            %     dimensionality j. Otherwise empty.
            %   loss
            %     Value of loss function for the training data. Computed
            %     using method specified in the params object.
            % Notes:
            %   The maximum likelihood PPCA covariance matrix is given in
            %   T&B section 3.2. Its eigenvectors are the eigenvectors of
            %   the sample covariance matrix. The first q eigenvalues are
            %   those of the sample covariance matrix. The remaining d-q
            %   eigenvalues are all sigma^2 (the maximum likelihood noise
            %   variance).
            
            % dimensionality of the observed data
            d = numel(mu);
            
            % number of nonzero eigenvalues
            tol = d * eps(eigval(1));
                % default tolerance as would be used by rank() on the
                % covariance matrix
            if eigval(end) > tol
                nnzeig = numel(eigval);
            else
                nnzeig = find(eigval > tol, 1, 'last');
            end
            
            % choose dimensionalty (q)
            
            % fixed dimensionality
            if obj.params.data.dim_criterion >= 1
                q = round(obj.params.data.dim_criterion);
                evidence = [];
            
            % choose based on fraction of variance accounted for
            elseif obj.params.data.dim_criterion > 0
                q = find( ...
                    cumsum(eigval) ./ sum(eigval) ...
                    > obj.params.data.dim_criterion, ...
                    1 ...
                );
                evidence = [];
            
            % use bayesian model selection
            else
                % compute log evidence for each dimensionality
                evidence = ...
                    obj.get_evidence(n, d, eigval(1:nnzeig), nnzeig-1);
                
                % choose dimensionality that maximizes it
                [max_le, q] = max(evidence);
            end
            
            % dimensionality must be less than the number of nonzero
            % eigenvalues
            q = min(q, nnzeig - 1);
            
            % maximum likelihood PPCA noise variance
            sigma2 = sum(eigval(q+1 : end)) ./ (d-q);
                % mean of the smallest/discarded eigenvalues
                % T&B equation 8
            
            % log determinant of the PPCA covariance matrix
            %eigval(1:q) = ones(size(eigval(1:q))) * mean(eigval(1:q));
            ldC = sum(log(eigval(1:q))) + (d-q)*log(sigma2);
                % determinant is the product of the eigenvalues
            
            % sum of first two terms for the log density of a data point
            ldt = -d/2*log(2*pi) - ldC/2;
            
            % compute value of loss function
            
            % if splitting based on bayesian evidence, return average
            % negative log evidence as the loss
            if obj.params.data.split_evidence
                if obj.params.data.dim_criterion >= 0
                    error([ ...
                        'Can''t split based on evidence unless ', ...
                        'also using it to choose dimensionality' ...
                    ]);
                end
                loss = -max_le ./ n;
            
            % otherwise loss is the average negative log likelihood of the
            % training set (i.e. cross entropy between the empirical
            % distribution and the model)
            else
                loss = d/2 - ldt;
                
                % for maximum likelihood PPCA model, the log likelihood of
                % the training simplifies to:
                %   -n*d/2*log(2*pi) - n/2*log(det(C)) - 1/2*n*d
                % where n, d are the number of points/dimensions, C is the
                % ppca covariance matrix
            end
        end % fit_ppca()
        
        
        % Compute Bayesian evidence for each choice of dimensionality for
        % PPCA model
        function [ evidence ] = get_evidence( obj, N, d, lambda, kmax )
            % Inputs:
            %   N
            %     Number of data points
            %   d
            %     Dimensionality of data
            %   lambda
            %     Eigenvalues of the sample covariance matrix, sorted in
            %     descending order. Zero eigenvalues must be discarded.
            %   kmax
            %     Maximum dimensionality to consider. Must be less than
            %     the number of nonzero eigenvalues.
            % Outputs:
            %   evidence
            %     evidence(j) contains the evidence for dimensionality j
            % Notes:
            %   Uses equations from:
            %   Minka (2000). Automatic choice of dimensionality for PCA.
            
            % pre-compute and cache reusable terms.
            % store in obj.params. will be accessible to all nodes that
            % share params object
            if ~all(isfield(obj.params.data, {'ev_m', 'ev_t125'}))
                
                % all possible choices of k
                all_k = 1:d-1;
                
                % parameter manifold dimensionality
                % (just above equation 52)
                all_m = d .* all_k - all_k .* (all_k + 1) ./ 2;
                
                % log p(U) (equation 52)
                lpU = ( ...
                    -all_k .* log(2) ...
                    - (all_m + all_k) .* log(pi)/2 ...
                    + cumsum(gammaln((d-all_k+1)./2)) ...
                ); % this looks different than the equation because it can
                   % be re-expressed/simplified, and some vars become
                   % equivalent to m.
                
                % terms for log evidence (from equation 80)
                term1 = 3 .* all_k ./ 2 .* log(2);
                term2 = lpU;
                term5 = (all_m+all_k) ./ 2 .* log(2*pi);
                
                % cache values
                obj.params.data.ev_m = all_m;
                obj.params.data.ev_t125 = term1 + term2 + term5;
            end
            
            % number of nonzero eigenvalues
            dn = numel(lambda);
            
            % values of k to try
            ks = 1 : kmax;
            
            % temporary variables used to speed up further computations
            lambda_inv = 1 ./ lambda;
                % inverse of eigenvalues
            cs_log_lambda = cumsum(log(lambda));
            cs_log_lambda = cs_log_lambda(1:kmax)';
                % cs_log_lambda(k) = sum_{i=1:k} log(lambda(i))
                % format as row vector, for k = 1:kmax
            
            % parameter manifold dimensionality
            m = obj.params.data.ev_m(1:kmax);
                % retrieve cached values
            
            % noise variance
            % (equation 81)
            tmp = cumsum(lambda(end : -1 : 1))';
            vhat = tmp(end-1 : -1 : 1) ./ (d - ks);
                % vhat(k) = 1/(d-k) sum_{i=k+1:d} lambda(i)
                % mean of discarded eigenvalues. same as maximum likelihood
                % noise variance in probabilistic pca paper. zero
                % eigenvalues were removed, so take the sum and divide by
                % the total number.
            
            % log determinant of the Hessian matrix (log |A_Z|)
            % (equation 73, using equations 70,81 for lambda hat)
            
            % the way we compute it here is based on taking the log of
            % equation 73 and rewriting it as:
            %   log |A_Z| (for a given choice of k) =
            %     sum_{i=1:k} sum_{j=i+1:k} log(1/lambda(j) - 1/lambda(i))
            %     + (d-k) * sum_{i=1:k} log(1/vhat(k) - 1/lambda(i))
            %     + sum_{i=1:k} sum_{j=i+1:dn} log(lambda(i) - lambda(j))
            %     + (d-dn) * sum_{i=1:k} log(lambda(i))
            %     + m*log(N)
            
            tmp1 = triu( ...
                log(bsxfun(@minus, lambda_inv(:)', lambda_inv(:))), ...
                1 ...
            );
            tmp1 = cumsum(sum(tmp1(1 : kmax, :), 1));
            tmp1 = tmp1(1:kmax);
                % tmp1(k) =
                % sum_{i=1:k} sum_{j=i+1:k} log(1/lambda(j) - 1/lambda(i))
            
            tmp2 = sum(tril(log(bsxfun(@minus, 1 ./ vhat(:), lambda_inv(:)'))), 2);
            tmp2 = (d - ks) .* tmp2';
                % tmp2(k) =
                % (d-k) * sum_{i=1:k} log(1/vhat(k) - 1/lambda(i));
            
            tmp3 = sum( ...
                triu(log(bsxfun(@minus, lambda(:), lambda(:)')), 1), ...
                2 ...
            );
            tmp3 = cumsum(tmp3);
            tmp3 = tmp3(1:kmax)';
                % tmp3(k) =
                % sum_{i=1:k} sum_{j=i+1:dn} log(lambda(i) - lambda(j))

            tmp4 = (d - dn) .* cs_log_lambda;
                % tmp4(k) = (d-dn) * sum_{i=1:k} log(lambda(i))

            tmp5 = m .* log(N);

            ldAZ = tmp1 + tmp2 + tmp3 + tmp4 + tmp5;
            
            % sum of terms 1, 2, 5 for log evidence
            t125 = obj.params.data.ev_t125(1:kmax);
                % retrieve cached values

            % log evidence (equation 80)
            term3 = -N/2 .* cs_log_lambda;
            term4 = -N .* (d-ks) ./ 2 .* log(vhat);
            term6 = -ldAZ ./ 2;
            term7 = -ks ./ 2 .* log(N);
            evidence = t125 + term3 + term4 + term6 + term7;
        end
        
        
        % Find optimal splitting threshold along a single direction
       function [ t, child_loss, loss ] = split_opt_1d( obj, x, y, nodes )
            % Inputs:
            %   x, y
            %     Sorted input/outputs. Rows correspond to points, columns
            %     to dimensions. The data is sorted in order of the
            %     projections of the inputs onto the normal vector of the
            %     splitting hyperplane.
            %   nodes
            %     A vector [node_L, node_R] containing node objects that
            %     will serve as the left/right child nodes produced by the
            %     split. Will be updated with trained local model
            %     parameters. No other object properties will be set.
            % Outputs:
            %   t
            %     Splitting threshold index. Sorted points 1:t are assigned
            %     to the left node and points t+1:end are assigned to the
            %     right node.
            %   child_loss
            %     A vector [loss_L, loss_R] giving the value of the loss
            %     function for the local model in each new child node on
            %     its assigned training points.
            %   loss
            %     Value of loss function for the computed split. Given by
            %     the average loss of the local model of the left/right
            %     child nodes, weighted by the number of training points
            %     assigned to each.
            % Notes:
            %   Attempts to bracket the optimum splitting threshold using
            %   an initial grid search. Then seeks the optimum in this
            %   interval using golden section search. The loss function may
            %   have multiple local minima. This strategy tries to find a
            %   good value without exhaustively searching every
            %   possibility.
            
            % number of points, dimensions
            [n, d] = size(y);
    
            % initial grid search
            
            % construct grid
            tmin = obj.params.data.min_leaf_pts;
            tmax = n - tmin;
                % no child node receives fewer than the minimum number of
                % points (this could be violated if many points share the
                % same projections)
            ngrid = min(max(obj.params.data.ngrid, 3), tmax-tmin+1);
            gt = unique(round(linspace(tmin, tmax, ngrid)));
                % each grid point is a potential splitting threshold index
            
            % best split found so far
            % initialize at first grid point
            [models, loss] = obj.update_split_models( ...
                y(1 : gt(1), :), y(gt(1)+1:end, :) ...
            ); % ppca models for data to left/right of split. rather than
            % recomputing from scratch for each threshold, update by
            % transfering points between left/right models.
            
            % value of loss function for each grid point
            gl = zeros(size(gt));
            gl(1) = loss;
            
            % loop thru grid points, compute loss for each
            my_models = models;
            for it = 2 : ngrid

                % add data to left model, remove from right model
                y_block = y(gt(it-1)+1 : gt(it), :);
                [my_models, my_loss] = ...
                    obj.update_split_models(my_models, y_block, 1);
                
                gl(it) = my_loss;
                
                % save models if they're the best so far
                if my_loss < loss
                    loss = my_loss;
                    models = my_models;
                end
            end % grid search
            
            % bracket the minimum
            % b is the grid point with minimum loss
            % a and c are grid points to the left/right of b
            % (or b +/- 1 if b is at the edge of the grid)
            [minval, minloc] = min(gl);
            b = gt(minloc);
            lb = minval;
            if minloc > 1
                a = gt(minloc - 1);
            else
                a = b - 1;
            end
            if minloc < ngrid
                c = gt(minloc + 1);
            else
                c = b + 1;
            end
            
            % a, b, c are splitting threshold indices
            % we now have a < b < c, loss(b) <= loss(a), loss(b) <= loss(c)
            % so interval [a, c] contains a local minimum of loss function
            % use golden section search to find it

            % the golden ratio
            w = (3-sqrt(5))/2;
            
            % DEBUG
            nevals = ngrid;

            % loop until no points remain betweeen a, b, c
            while c - a > 2

                % choose a new point d in the interval (a, b) or (b, c)
                % (whichever is wider)
                if (c-b) >= (b-a)
                    d = round(b + w*(c-b));
                else
                    d = round(b - w*(b-a));
                end
                
                % data to be transfered between left/right models
                pt1 = min(b, d) + 1;
                pt2 = max(b, d);
                y_block = y(pt1 : pt2, :);
                
                % update models, evaluate loss function at new point
                [my_models, ld] = ...
                    obj.update_split_models(models, y_block, d - b);
                    % if d > b, add data to left model, remove from right
                    % otherise do opposite
                
                % DEBUG
                nevals = nevals + 1;

                % three of the four points will become the triplet for the
                % next iteration, respecting the constraints that
                % a < b < c, loss(b) <= loss(a), loss(b) <= loss(c)
                % models struct always corresponds to point b
                if d > b
                    if ld >= lb
                        c = d;
                    else
                        a = b;
                        b = d;
                        lb = ld;
                        models = my_models;
                    end
                else
                    if ld >= lb
                        a = d;
                    else
                        c = b;
                        b = d;
                        lb = ld;
                        models = my_models;
                    end
                end
            end % golden section search

            % outputs
            
            % optimal splitting threshold index, value of loss function
            t = b;
            loss = lb;
            
            % value of loss function for left/right ppca model
            child_loss = [models.loss_L, models.loss_R];
            
            % update nodes with trained local models
            q_L = models.q_L;
            nodes(1).d = obj.d;
            nodes(1).q = q_L;
            nodes(1).mu = models.spca.mu_L;
            nodes(1).U = models.spca.eigvec_L(:, 1:q_L);
            nodes(1).lambda = models.spca.eigval_L(1:q_L);
            nodes(1).sigma2 = models.sigma2_L;
            nodes(1).ldt = models.ldt_L;
            nodes(1).eigs = models.spca.eigval_L;
            nodes(1).evidence = models.evidence_L;
            
            q_R = models.q_R;
            nodes(2).d = obj.d;
            nodes(2).q = q_R;
            nodes(2).mu = models.spca.mu_R;
            nodes(2).U = models.spca.eigvec_R(:, 1:q_R);
            nodes(2).lambda = models.spca.eigval_R(1:q_R);
            nodes(2).sigma2 = models.sigma2_R;
            nodes(2).ldt = models.ldt_R;
            nodes(2).eigs = models.spca.eigval_R;
            nodes(2).evidence = models.evidence_R;
        end % split_opt_1d()
        
        
        % Build and update PPCA models for a set of points partitioned into
        % left/right parts.
        function [ models, loss ] = update_split_models( obj, varargin )
            % Usage:
            %   update_split_models( obj, y_L, y_R )
            %     Used to initialize models
            %   update_split_models( obj, models, y, direction )
            %     Used to update models
            % Inputs:
            %   y_L, y_R
            %     Output points assigned to left/right models. Rows
            %     correspond to points, columns to dimensions.
            %   models
            %     Struct specifying models to be updated. See outputs for
            %     description.
            %   y
            %     Output points to transfer between left/right models. Rows
            %     correspond to points, columns to dimensions. Must be part
            %     of data set used to initialize models.
            %   direction
            %     A scalar. If >=0, points will be transfered from the
            %     right model to the left model. If negative, points will
            %     be transfered from left to right.
            % Outputs:
            %   models
            %     Struct containing parameters for left/right PPCA models.
            %     Contains fields:
            %       spca
            %         SplitPCA object. Contains number of points, means,
            %         eigendecomposition of the sample covariance matrix
            %         for left/right models.
            %       Other fields correspond to the outputs of fit_ppca()
            %   loss
            %     Average loss for the PPCA model in the left/right
            %     child nodes (weighted by the number of points in each).
            %     The loss function to use is specified in obj.params.
            % Notes:
            %   This function gives a way to update the models by
            %   transfering points between the left/right models, instead
            %   of recomputing from scratch.
            
            % initialize models
            if nargin == 3
                spca = SplitPCA().train(varargin{1}, varargin{2});
            
            % update models with new data block
            elseif nargin == 4
                
                models = varargin{1};
                y = varargin{2};
                direction = varargin{3};
                
                spca = models.spca.update(y, direction);
            else
                error('Usage error');
            end
            
            % fit ppca models
            [q_L, sigma2_L, ldt_L, evidence_L, loss_L] = obj.fit_ppca( ...
                spca.n_L, ...
                spca.mu_L, ...
                spca.eigval_L, ...
                spca.eigvec_L ...
            );
            [q_R, sigma2_R, ldt_R, evidence_R, loss_R] = obj.fit_ppca( ...
                spca.n_R, ...
                spca.mu_R, ...
                spca.eigval_R, ...
                spca.eigvec_R ...
            );
            
            % pack results into struct
            models = struct( ...
                'spca', spca, ...
                'q_L', q_L, 'q_R', q_R, ...
                'sigma2_L', sigma2_L, 'sigma2_R', sigma2_R, ...
                'ldt_L', ldt_L, 'ldt_R', ldt_R, ...
                'evidence_L', evidence_L, 'evidence_R', evidence_R, ...
                'loss_L', loss_L, 'loss_R', loss_R ...
            );
            
            % loss function
            nL = spca.n_L;
            nR = spca.n_R;
            loss = nL/(nL+nR)*loss_L + nR/(nL+nR)*loss_R;
        end % update_split_models()
        
        
        % Predict output, given input
        function [ yhat ] = predict( obj, x )
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions.
            % Outputs:
            %   yhat
            %     Matrix containing predicted output for each input. Rows
            %     correspond to points, columns to dimensions. Predicted
            %     output is the mean of the local PPCA model.
            % Preconditions:
            %   Node's local model must be trained.
            
            yhat = repmat(obj.mu, size(x, 1), 1);
        end % predict()
        
        
        % Evaluate PPCA model to get log density of a set of outputs
        function [ lp ] = log_density( obj, y )
            % Inputs:
            %  y
            %    A set of outputs. Rows corresponds to points, columns to
            %    dimensions.
            % Outputs:
            %   lp
            %     lp(i) gives the log density of the ith point
            
            if size(y, 2) ~= obj.d
                error('Data dimensionality doesn''t match model');
            end
            
            % center data, project onto top q eigenvectors
            yc = bsxfun(@minus, y, obj.mu);
            v = yc * obj.U;
            
            % expression for (y-mu)^T * inv(C) * (y-mu)
            tmp = sum(v * diag(1 ./ obj.lambda) .* v, 2);
            tmp = tmp + ( ...
                sum(yc .^ 2, 2) ...
                - sum(v .^ 2, 2) ...
            ) ./ obj.sigma2;
        
            % log density
            lp = obj.ldt - tmp ./ 2;
            
            % this corresponds to the expression:
            %   log(p(y)) =
            %     -d/2 * log(2*pi)
            %     - 1/2 * log(det(C))
            %     - 1/2 * (y-mu)^T * inv(C) * (y-mu)
            % where mu is the mean, C is the covariance matrix, and y is
            % a new point (represented as a column vector)
        end % log_density()
        
        
        % Evaluate PPCA model to get log density of a set of outputs
        function [ lp ] = lpygx( obj, x, y )
            % Inputs:
            %  y
            %    A set of outputs. Rows corresponds to points, columns to
            %    dimensions.
            % Outputs:
            %   lp
            %     lp(i) gives the log density of the ith point
            
            if size(y, 2) ~= obj.d
                error('Data dimensionality doesn''t match model');
            end
            
            % center data, project onto top q eigenvectors
            yc = bsxfun(@minus, y, obj.mu);

            v = yc * obj.U;
            
            % expression for (y-mu)^T * inv(C) * (y-mu)
            tmp = sum(v * diag(1 ./ obj.lambda) .* v, 2);
            tmp = tmp + ( ...
                sum(yc .^ 2, 2) ...
                - sum(v .^ 2, 2) ...
            ) ./ obj.sigma2;
        
            % log density
            lp = obj.ldt - tmp ./ 2;
            
            % this corresponds to the expression:
            %   log(p(y)) =
            %     -d/2 * log(2*pi)
            %     - 1/2 * log(det(C))
            %     - 1/2 * (y-mu)^T * inv(C) * (y-mu)
            % where mu is the mean, C is the covariance matrix, and y is
            % a new point (represented as a column vector)
        end % lpygx()
        
        
        % Get PPCA covariance matrix
        function [ C ] = cov( obj )
            % Outputs:
            %   C
            %     Maximum likelihood PPCA covariance matrix
            
            % corresponds to replacing bottom d-q eigenvalues of the
            % sample covariance matrix with their mean (i.e. maximum
            % likelihood noise variance)
            
            C = ( ...
                obj.U * diag(obj.lambda) * obj.U' ...
                + obj.sigma2 .* (eye(obj.d) - obj.U * obj.U') ...
            );
        end % cov()
    end % methods
end % classdef