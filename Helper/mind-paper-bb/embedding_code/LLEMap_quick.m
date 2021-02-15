classdef LLEMap_quick < handle
    
    properties
        Xtrain; % (matrix)
        % Input training points (points x dimensions)
        Ytrain; % (matrix)
        % Output training points (points x dimensions)
        h; % (real)
        % Kernel bandwidth, in units of distance
        h_scaled;
        % Kernel bandwidth
        all_h; % (vector)
        % All possible choices of h, scaled
        searchmode;
    end % properties
    
    methods
        
        % Constructor
        function [ obj ] = LLEMap_quick( h )
            % Inputs:
            %   k
            %     Number of nearest neighbors. If a vector of values is
            %     given, one will be chosen by cross validation.
            %   h (Optional)
            %     Regularization parameter. If a vector of values is given,
            %     one will be chosen by cross validation. [Default: 0]
            %   nfolds (Optional)
            %     Number of cross validation folds. [Default: 10]
            obj.all_h = h(:)';
            
        end % constructor
        
        
        % Fit mapping
        function [] = fit( obj, X, Y, searchmode)
            % Inputs:
            %   X
            %     Matrix of training input points (points x dimensions)
            %   Y
            %     Matrix of training output points (points x dimensions)
            % Posteffects:
            %   Stores X and Y as object properties. Chooses k and h
            %   using cross validation if multiple values were provided.
            
            if ~exist('searchmode','var')
                obj.searchmode = 'opt';
            else
                obj.searchmode = searchmode;
            end
            
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
            
            % set k and h if no cross validation requested
            nh = numel(obj.all_h);
            if nh == 1
                obj.h = obj.all_h;
                return;
            end
            
            % otherwise choose k and h by exact analytic loo cv
            D = squareform(pdist(X));
            all_h_scaled = obj.all_h * mean(D(:));
            
            % best values so far
            best_h = nan;
            best_err = inf;
            
            switch obj.searchmode
                case 'grid'
                    % loop thru h
                    for it_h = 1 : nh
                        my_h = all_h_scaled(it_h);
                        err = loo_error(my_h,D,Y);
                        % keep parameters if they improve cost
                        if err < best_err
                            best_h = my_h;
                            best_err = err;
                        end
                    end % loop thru h
                case 'opt'
                    hbnds = [min(all_h_scaled) max(all_h_scaled)];
                    options = optimset('TolX',10^-3,'MaxIter',50);
                    errfun = @(my_h) loo_error(my_h,D,Y);
                    best_h = fminbnd(errfun,hbnds(1),hbnds(2),options);
            end
            % store best parameters
            obj.h = best_h;
            obj.h_scaled = obj.h / mean(D(:));
        end % fit()
        
        % Map new input points into output space
        function [ Y ] = transform( obj, X )
            
            % check inputs
            if isempty(obj.Xtrain)
                error('Model must be fit first');
            end
            if size(X, 2) ~= size(obj.Xtrain, 2)
                error('Dimensionality of points doesn''t match model');
            end
            
            % compute weight matrix
            D = pdist2(obj.Xtrain,X);
            L = exp(-(D/obj.h).^2);
            L = L ./ sum(L,1);
            
            % compute fit
            Y = (obj.Ytrain' * L)';
        end % transform()
        
    end % methods
end % classdef

function err = loo_error(my_h, D, Y)
    L = exp(-(D/my_h).^2);
    L = L ./ sum(L,1);
    Yhat = (Y' * L)';
    err = mean((norm_rows(Yhat - Y) ./ (1 - diag(L))).^2);
end


