classdef SplitPCA
    properties
        d; % (scalar)
            % Dimensionality of input data
        n_L; % (scalar)
        n_R;
            % Number of data points assigned to left/right model
        s_L; % (1 x d) vector
        s_R;
            % Sum of data points for left/right model
        O_L; % (d x d matrix)
        O_R;
            % Sum of outer products of data points for left/right model
        mu_L; % (1 x d vector)
        mu_R;
            % Sample mean for left/right model
        C_L; % (d x d matrix)
        C_R;
            % Sample covariance matrix for left/right model
        eigval_L; % (r x 1 vector)
        eigval_R;
            % Eigenvalues of the covariance matrix for left/right model
            % Sorted in descending order
            % r is the rank of the covariance matrix
        eigvec_L; % (d x r matrix)
        eigvec_R;
            % Corresponding eigenvectors of the covariance matrix, stored
            % on the columns. For left/right model.
    end % properties
    
    methods
        function [ obj ] = train( obj, x_L, x_R )
            
            % number of points, dimensions
            [obj.n_L, obj.d] = size(x_L);
            obj.n_R = size(x_R, 1);
            
            % sum of data points
            obj.s_L = ones(1, obj.n_L) * x_L;
            obj.s_R = ones(1, obj.n_R) * x_R;
            
            % sum of outer products
            obj.O_L = x_L' * x_L;
            obj.O_R = x_R' * x_R;
            
            % compute mean, covariance matrix, eigendecomposition
            obj = obj.compute_eig();
        end % train()
        
        
        function [ obj ] = update( obj, x, direction )
            
            % number of points, dimensions
            [n, d] = size(x);
            
            % sum of data points, sum of outer products
            s = ones(1, n) * x;
            O = x' * x;
            
            % add points to left model, remove from right model
            if direction >= 0
                obj.n_L = obj.n_L + n;
                obj.n_R = obj.n_R - n;
                obj.s_L = obj.s_L + s;
                obj.s_R = obj.s_R - s;
                obj.O_L = obj.O_L + O;
                obj.O_R = obj.O_R - O;
            % remove points from left model, add to right model
            else
                obj.n_L = obj.n_L - n;
                obj.n_R = obj.n_R + n;
                obj.s_L = obj.s_L - s;
                obj.s_R = obj.s_R + s;
                obj.O_L = obj.O_L - O;
                obj.O_R = obj.O_R + O;
            end
            
            % compute mean, covariance matrix, eigendecomposition
            obj = obj.compute_eig();
        end % update()
        
        
        function [ obj ] = compute_eig( obj )
            
            % sample mean
            obj.mu_L = obj.s_L ./ obj.n_L;
            obj.mu_R = obj.s_R ./ obj.n_R;
            
            % sample covariance matrix
            obj.C_L = obj.O_L ./ obj.n_L - obj.mu_L' * obj.mu_L;
            obj.C_R = obj.O_R ./ obj.n_R - obj.mu_R' * obj.mu_R;
            
            % eigendecomposition of covariance matrix
            [eigvec, eigval] = eig(obj.C_L, 'vector');
            [obj.eigval_L, ord] = sort(eigval, 'descend');
            obj.eigvec_L = eigvec(:, ord);
            [eigvec, eigval] = eig(obj.C_R, 'vector');
            [obj.eigval_R, ord] = sort(eigval, 'descend');
            obj.eigvec_R = eigvec(:, ord);
        end % compute_eig()
    end % methods
end % classdef