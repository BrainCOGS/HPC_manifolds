% DForest class
% Implements decision forests composed of DTrees

classdef DForest < handle
    
    properties
        
        % Forest parameters
        ntrees; % (integer)
            % How many trees in the forest
        verbose; % (logical)
            % Whether to print progress when training
        callback; % (function handle)
            % Function to execute after each tree is trained.
            % Must must have the form: f(obj, x, y, n)
            % Where obj is the trained tree object, x and y are the data
            % used to train the tree, and n is the index of the tree in
            % the forest.
        
        % Tree parameters
        patience;
            % How many times we're willing to attempt splitting a leaf node
            % before giving up on it. An attempt succeeds if it reduces
            % the overall loss.
        
        % Model
        trees; % (array of DTree)
            % All trees in the forest
        ndims_in; % (integer)
            % Dimensionality of input
    end
    
    
    methods
        
        % Constructor
        function [ obj ] = DForest( node, varargin )
            % Inputs:
            %   node (TreeNode object)
            %     Prototype for the root node of all trees. Determines what
            %     type of decision tree to create (classification,
            %     regrssion, etc.). Can be obtained by calling the
            %     constructor of the the appropriate TreeNode subclass. All
            %     root nodes will be created as objects of the same class,
            %     and node's shared parameters will be copied over.
            %   Further inputs given as optional name/value pairs. See
            %   class properties for names, descriptions. Default values
            %   in constructor code.
            % Notes:
            %   All nodes in all trees/forests constructed from the same
            %   prototype node will share the same underlying 'params'
            %   object.
            
            % default parameters
            default_args = struct( ...
                'ntrees', 1, ...
                'verbose', false, ...
                'callback', [], ...
                'patience', 1 ...
            );
            
            % update with user supplied args, set object properties
            args = parse_args(varargin, default_args, false);
            obj.ntrees = args.ntrees;
            obj.verbose = args.verbose;
            obj.callback = args.callback;
            obj.patience = args.patience;
            
            % initialize trees
            obj.trees = repmat(DTree(node), 1, obj.ntrees);
            for it = 1 : obj.ntrees
                obj.trees(it) = DTree( ...
                    node, ...
                    'patience', obj.patience ...
                );
            end
        end % constructor
        
        
        % Train forest
        function [] = train( obj, x, y )
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions
            %   y
            %     Output data matrix. Rows correspond to points, columns to
            %     dimensions
            % Postconditions:
            %   All trees in forest are trained.
            
            % number of dimenions/points
            if ~ismatrix(x) || ~ismatrix(y) || size(x, 1) ~= size(y, 1)
                error('Invalid training data');
            end
            obj.ndims_in = size(x, 2);
            
            % train each tree
            %parfor it = 1 : obj.ntrees
            for it = 1 : obj.ntrees
                
                % print progress if requested
                if obj.verbose
                    fprintf('%d/%d\n', it, obj.ntrees);
                end
                
                % train current tree
                obj.trees(it).train(x, y);
                
                % execute callback function
                if ~isempty(obj.callback)
                    obj.callback(obj.trees(it), x, y, it);
                end
            end
        end
        
        
        % Get predicted output for given inputs
        function [ yhat, all_ids, all_yhat ] = predict( obj, x, max_depth )
            % Inputs:
            %   x
            %     Data matrix. Rows correspond to points, columns to
            %     dimensions.
            %   max_depth (Optional)
            %     Maximum depth to traverse in each tree before before
            %     assigning each data point to a node and using its local
            %     model to compute predicted output. If given as a vector,
            %     results will be returned for each value (must be sorted
            %     in ascending order). [Default: Traverse all the way to
            %     leaf nodes].
            % Outputs:
            %   yhat
            %     Predicted output. Rows correspond to points, columns to
            %     dimensions. Third dimension corresponds to different
            %     choices of max_depth, if multiple values given.
            %   ids
            %     id of node to which each data point is assigned. Rows
            %     correspond to points. Columns correspond to different
            %     choices of max_depth, if multiple values given.
            % Notes:
            %   The max_depth parameter allows a deep tree to be trained,
            %   but evaluated at shallower depths. Passing multiple values
            %   allows for efficiently testing multiple depths, because
            %   only a single traversal of the tree is required.
            
            % check that model has been trained
            if ( ...
                    any([obj.trees.nnodes] == 0) ...
                    || isempty(obj.trees(end).root.id) ...
                )
                error('Model must be trained first');
            end
            
            % input dimensions
            [npts, ndims_in] = size(x);
            if ndims_in ~= obj.ndims_in
                error('Input dimensions don''t match model');
            end
            
            % predict using leaf nodes by default
            if nargin < 3 || isempty(max_depth)
                max_depth = Inf;
            end
            
            % how many values of max_depth to try
            ndepths = numel(max_depth);
            if ~issorted(max_depth)
                error('Values in max_depth must be in ascending order');
            end
            if max_depth(1) < 0
                error('Invalid max_depth');
            end
            
            % predict first data point using first tree's root node to
            % determine the number of output dimensions
            tmp = obj.trees(1).root.predict(x(1, :));
            ndims_out = size(tmp, 2);
            
            % initialize outputs
            yhat = zeros(npts, ndims_out, ndepths);
            if nargout >= 2
                all_ids = zeros(npts, obj.ntrees, ndepths);
            end
            if nargout >= 3
                all_yhat = zeros(npts, ndims_out, obj.ntrees, ndepths);
            end
            
            % loop thru trees
            for it = 1 : obj.ntrees
                [my_yhat, my_ids] = obj.trees(it).predict(x, max_depth);
                yhat = yhat + my_yhat;
                if nargout >= 2
                    all_ids(:, it, :) = my_ids;
                end
                if nargout >= 3
                    all_yhat(:, :, it, :) = my_yhat;
                end
            end
            yhat = yhat ./ obj.ntrees;
        end % predict()
        
    end % methods
end % classdef