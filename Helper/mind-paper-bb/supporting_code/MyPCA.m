% MyPCA()
% Principal component analysis

classdef MyPCA < handle
    
    properties
        mu;
            % Mean vector
        W;
            % Weight matrix (input dimensions x components)
            % (Eigenvectors of the covariance matrix)
        lambda;
            % Eigenvalues of the covariance matrix
            % (Variance of projections onto each weight vector)
    end
    
    methods
        
        % Train model
        function [ Y ] = fit( obj, X, k )
            % Inputs:
            %   X
            %     Training data (points x dimensions)
            %   k (optional)
            %     Number of components to include in transformed data. If
            %     0<k<1, include top components that explain at least this
            %     fraction of the variance. [Default: all]
            % Outputs:
            %   Y
            %     Transformed data
            
            % number of points, dimensions
            [n, d] = size(X);
            
            % center data
            obj.mu = mean(X, 1);
            X = bsxfun(@minus, X, obj.mu);
            
            % find weights, variance of projections
            % choose faster method depending on number of points vs. dims
            
            % eigendecomposition of covariance matrix
            if n > d
                [V, lambda] = eig(cov(X, 1), 'vector');
                [lambda, idx] = sort(lambda, 'descend');
                V = V(:, idx);
                s = sqrt(n*lambda);
            
            % svd of data matrix
            else
                [~, S, V] = svd(X, 'econ');
                s = diag(S);
                lambda = s.^2 / n;
            end
            
            % exclude components with zero variance
            tol = max(n, d) * eps(max(s));
            sel = s > tol;
                % same method used by rank()
            
            % store weights, variances
            obj.W = V(:, sel);
            obj.lambda = lambda(sel);
            
            % return transformed data if requested
            if nargout > 0
                if nargin < 3 || isempty(k)
                    k = numel(lambda);
                elseif k > 0 && k < 1
                    k = find(obj.variance_explained() >= k, 1);
                else
                    k = min(round(k), numel(lambda));
                end
                Y = X * V(:, 1:k);
            end
        end % fit()
        
        
        % Reduce dimensionality of high dimensional points
        function [ Y ] = transform( obj, X, k )
            % Inputs:
            %   X
            %     Data (points x dimensions)
            %   k (optional)
            %     Number of components to include in transformed data. If
            %     0<k<1, include top components that explain at least this
            %     fraction of the variance. [Default: all]
            % Outputs:
            %   Y
            %     Transformed data (points x components). Computed by
            %     subtracting the (model) mean, then projecting onto the
            %     weight vectors of the top k components.
            
            % check input
            if isempty(obj.mu)
                error('Model has not been fit');
            end
            if size(X, 2) ~= numel(obj.mu)
                error('Dimensionality of data doesn''t match model');
            end
            
            % choose number of components
            if nargin < 3 || isempty(k)
                k = numel(obj.lambda);
            elseif k > 0 && k < 1
                k = find(obj.variance_explained() >= k, 1);
            else
                k = min(round(k), numel(obj.lambda));
            end
            
            % transform data
            Y = bsxfun(@minus, X, obj.mu) * obj.W(:, 1:k);
        end % transform()
        
        
        % Compute fraction of variance explained
        function [ r2 ] = variance_explained( obj )
            % Outputs:
            %   r2
            %     r2(k) gives the fraction of variance explained (R^2) by
            %     the top k components
            
            if isempty(obj.lambda)
                error('Model has not been fit');
            end
            
            r2 = cumsum(obj.lambda) ./ sum(obj.lambda);
        end % variance_explained()
        
        
        % Reconstruct high dimensional points from low dimensional
        % projections
        function [ X ] = inverse_transform( obj, Y )
            % Inputs:
            %   Y
            %     Data (n points x k components). Represents projections
            %     along the top k principal components.
            % Outputs:
            %   X
            %     Reconstructed points in input space (points x dimensions)
            
            % number of points, components
            [n, k] = size(Y);
            
            % check input
            if isempty(obj.mu)
                error('Model has not been fit');
            end
            if k > numel(obj.lambda)
                error('Dimensionality of data doesn''t match model');
            end
            
            % perform inverse transform
            X = bsxfun(@plus, Y * obj.W(:, 1:k)', obj.mu);
        end

    end % methods
end % classdef