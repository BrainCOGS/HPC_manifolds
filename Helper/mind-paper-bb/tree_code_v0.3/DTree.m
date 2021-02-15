% DTree class
% Decision tree class. Can implement classification trees, regression
% trees, etc.

classdef DTree < handle
    
    properties
        % Parameters
        patience;
            % How many times we're willing to attempt splitting a leaf node
            % before giving up on it. An attempt succeeds if it reduces
            % the overall loss.
        
        % Model
        root; % (TreeNode object)
            % The root node
        ndims_in; % (integer)
            % Dimensionality of input
        
        % Tree stats
        nnodes = 0; % (integer)
            % How many nodes in tree
        nleaves = 0; % (integer)
            % How many leaves in tree
        depth = 0; % (integer)
            % Maximum depth of any leaf node (root node has depth 0)
    
    end % properties
    
    methods
        
        % Constructor
        function [ obj ] = DTree( node, varargin )
            % Inputs:
            %   node (TreeNode object)
            %     Prototype for the tree's root node. Determines what type
            %     of decision tree to create (classification, regression,
            %     etc.). Obtained by calling the constructor of the the
            %     appropriate TreeNode subclass. The root will be created
            %     as an object of the same class, and node's shared
            %     parameters will be copied to the root.
            %   Further inputs given as optional name/value pairs:
            %     'patience'
            %       Maximum number of times to try splitting a given node.
            %       node. If greater than 1, multiple attempts will be made
            %       if the first fails (e.g. doesn't reduce the loss
            %       function). [Default: 1]
            % Notes:
            %   If multiple trees are created from the same 'node' input
            %   object, all nodes across all trees will share the same
            %   underlying 'params' object.
            
            if nargin == 0
                return;
            end
            
            % create root node
            if isa(node, 'TreeNode')
                obj.root = node.spawn();
                    % get a new node of the same type as prototype.
                    % copies node's shared params object to root.
                obj.root.depth = 0;
            else
                error('node must be a TreeNode object');
            end
            
            % default parameters
            default_args = struct( ...
                'patience', 1 ...
            );
            
            % update with user supplied args, set object properties
            args = parse_args(varargin, default_args, false);
            obj.patience = args.patience;
        end % constructor
        
        
        % Train tree
        function [] = train( obj, x, y )
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions
            %   y
            %     Output data matrix. Rows correspond to points, columns to
            %     dimensions
            % Postconditions:
            %   Entire tree is trained, consisting of new nodes linked to
            %   root node. All nodes have a trained local model. All non-
            %   leaf nodes have splitting parameters set.
            % Notes:
            %   If the tree is alreaady trained, the current model will be
            %   destroyed and a new one trained.
            
            % number of dimenions/points
            if ~ismatrix(x) || ~ismatrix(y) || size(x, 1) ~= size(y, 1)
                error('Invalid training data');
            end
            obj.ndims_in = size(x, 2);
            
            % delete root node's connections if they exist
            % (i.e. start from an empty tree)
            obj.root.parent = [];
            obj.root.children = [];
            obj.root.id = [];
            
            % train root node's local model
            obj.root.train(x, y);
            
            % set number of nodes/leaves, depth of tree
            obj.nnodes = 1;
            obj.nleaves = 1;
            obj.depth = 0;
            
            % recursively split/train nodes
            obj.recursive_split(obj.root, x, y);
            
            % assign id to every node in tree
            obj.label_nodes();
        end % train()
        
        
        % Recursively split/train nodes
        function [] = recursive_split( obj, node, x, y )
            % Inputs:
            %   node
            %     Current leaf node to split (TreeNode object). Must
            %     already contain a trained local model and no children.
            %   x, y
            %     Input/output training points assigned to current
            %     node. Rows correspond to points, columns to
            %     dimensions.
            % Notes / postconditions:
            %   Attempts to split node and train local models in the
            %   resulting child nodes. Recursively repeats this process in
            %   a depth-first manner until termination conditions are met.
            
            % stop if the current node isn't splittable
            % (too few points, or already at max depth or min loss)
            if ~node.is_splittable()
                return;
            end
            
            % try to split current node
            % if requested, try multiple times if unsuccessful
            success = false;
            attempts_left = obj.patience;
            while attempts_left > 0
                as = node.split(x, y);
                if ~isempty(as)
                    success = true;
                    attempts_left = 0;
                else
                    attempts_left = attempts_left - 1;
                end
            end
            
            % if split successful
            if success
                % update tree stats
                nchildren = size(as, 2);
                obj.nnodes = obj.nnodes + nchildren;
                obj.nleaves = obj.nleaves + nchildren - 1;
                obj.depth = max(obj.depth, node.depth + 1);
                
                % recursively split/train children
                for it = 1 : nchildren
                    obj.recursive_split( ...
                        node.children(it), ...
                        x(as(:, it), :), ...
                        y(as(:, it), :) ...
                    );
                end
            end
        end % recursive_split()
        
        
        % Assign id to every node in the tree
        function [] = label_nodes( obj )
            % Postconditions:
            %   The 'id' property of every node in the tree will be set to
            %   a unique integer. Ids are assigned incrementally in a
            %   breadth-first manner, starting with 1 at the root.
            
            % construct a list that will eventually contain all nodes
            q = obj.root; % add root node as the first item
            q(obj.nnodes) = obj.root;
                % only used to create an array w/ proper type
            last = 1; % last assigned item in the list
            
            % loop thru nodes
            for it = 1 : obj.nnodes
                % assign id to next node in list
                node = q(it);
                node.id = it;
                
                % append current node's children to end of list
                nchildren = numel(node.children);
                if nchildren > 0
                    q(last + 1 : last + nchildren) = node.children;
                    last = last + nchildren;
                end
            end
        end
        
        
        % Get predicted output for given inputs
        function [ yhat, ids ] = predict( obj, x, max_depth )
            % Inputs:
            %   x
            %     Data matrix. Rows correspond to points, columns to
            %     dimensions.
            %   max_depth (Optional)
            %     Maximum depth to traverse in tree before before assigning
            %     each data point to a node and using its local model to
            %     compute predicted output. If given as a vector, results
            %     will be returned for each value (must be sorted in
            %     ascending order). [Default: Traverse all the way to
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
            if obj.nnodes == 0 || isempty(obj.root.id)
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
            
            % predict first data point using root node to determine the
            % number of output dimensions
            tmp = obj.root.predict(x(1, :));
            ndims_out = size(tmp, 2);
            
            % outputs
            yhat = nan(npts, ndims_out, ndepths);
            ids = zeros(npts, ndepths);
            
            % recursively traverse tree to assign points to nodes and
            % compute predictions
            do_predict(obj.root, 1 : npts, 1);
            
            
            function [] = do_predict( node, idx, mdi )
                % Inputs:
                %   node
                %     Current node (TreeNode object)
                %   idx
                %     Indices of data points in current node
                %   mdi
                %     Index into max_depth vector, specifying next max
                %     depth parameter to work with. Model will be evaluated
                %     if depth is just about to exceed max_depth(mdi), or
                %     if a leaf node is hit.
                % Postconditions:
                %   yhat and ids will be updated.
                % Notes
                %   Called recursively. Terminates if leaf node or final
                %   max depth is hit.
                
                % is current node a leaf?
                is_leaf = isempty(node.children);
                
                % if this is a leaf node
                % or we're just about to exceed next max_depth
                if is_leaf || node.depth + 1 > max_depth(mdi)
                    
                    % use local model to predict output for current inputs
                    my_yhat = node.predict(x(idx, :));
                    
                    % if this is a leaf node, assign output to all further
                    % values of max_depth
                    if is_leaf
                        for it = mdi : ndepths
                            yhat(idx, :, it) = my_yhat;
                            ids(idx, it) = node.id;
                        end
                    
                    % otherwise assign it to current max_depth
                    else
                        yhat(idx, :, mdi) = my_yhat;
                        ids(idx, mdi) = node.id;
                    end
                    
                    % move on to next value in max_depth
                    mdi = mdi + 1;
                end
                
                % if this is a non-leaf node and we haven't exceeded final
                % maximum depth
                if ~is_leaf && mdi <= ndepths

                    % assign points to child nodes
                    as = node.assign(x(idx, :));
                    
                    % recursively pass data to children for prediction
                    nchildren = size(as, 2);
                    for it = 1 : nchildren
                        do_predict( ...
                            node.children(it), ...
                            idx(as(:, it)), ...
                            mdi ...
                        );
                    end
                end
            end % do_predict()
        end % predict()
        
        
        % Compute the log probability of outputs, given inputs
        function [ lpy ] = lpygx( obj, x, y, all )
            
            % check that node supports lpygx() method
            if ~ismethod(obj.root, 'lpygx')
                error(sprintf( ...
                    'Node type %s doesn''t support this operation', ...
                    class(obj.root) ...
                ));
            end
            
            % by default, compute log probability of each point in y given
            % corresponding point in x
            if nargin < 4
                all = false;
            end
            
            % number of points, dimensions
            [npts_x, ndims_x] = size(x);
            [npts_y, ndims_y] = size(y);
            if ndims_x ~= obj.ndims_in
                error('Input dimensions don''t match model');
            end
            
            % initialize output
            if all
                lpy = zeros(npts_x, npts_y);
            else
                if npts_x ~= npts_y
                    error([ ...
                        'x and y must contain the same ', ...
                        'number of points' ...
                    ]);
                end
                lpy = zeros(npts_x, 1);
            end
            
            doit(obj.root, 1 : npts_x);
            
            function [] = doit( node, idx )
                if isempty(idx)
                    % do nothing
                elseif isempty(node.children)
                    if all
                        for it = idx(:)'
                            lpy(it, :) = node.lpygx(x(idx, :), y )';
                        end
                    else
                        lpy(idx) = node.lpygx(x(idx, :), y(idx, :) );
                    end
                else
                    % assign points to child nodes
                    as = node.assign(x(idx, :));
                    
                    % recursively pass data to children
                    nchildren = size(as, 2);
                    for it = 1 : nchildren
                        doit(node.children(it), idx(as(:, it)));
                    end
                end
            end
        end
        

        % Get nodes containing given input points
        function [ nodes ] = get_nodes( obj, x )
            
            % input dimensions
            [npts, ndims_in] = size(x);
            
            % initialize node list
            nodes = repmat(obj.root.spawn(), npts, 1);
            
            % recursively traverse tree to assign points to nodes
            doit(obj.root, 1 : npts);
            
            function [] = doit( node, idx )
                if isempty(node.children)
                    nodes(idx) = node;
                else
                    % assign points to child nodes
                    as = node.assign(x(idx, :));
                    
                    % recursively pass data to children for prediction
                    nchildren = size(as, 2);
                    for it = 1 : nchildren
                        doit(node.children(it), idx(as(:, it)));
                    end
                end
            end

        end % get_nodes()
        
        
        % Get vector containing all leaf nodes
        function [ leaves ] = get_all_leaves( obj )
            
            % initialize node list
            leaves = repmat(obj.root.spawn(), obj.nleaves, 1);
            ptr = 1;
            
            % recursively traverse tree to assign points to nodes
            doit(obj.root);
            
            function [] = doit( node )
                children = node.children;
                
                % add node to list if it's a leaf
                if isempty(children)
                    leaves(ptr) = node;
                    ptr = ptr + 1;
                
                % otherwise, recursively apply function to children
                else
                    for it = 1 : numel(children)
                        doit(children(it));
                    end
                end
            end

        end % get_all_leaves()
        
        
        % Plot tree
        function [] = plot( obj, cvalfun )
            % Plots graph of the tree object. X coordinates are in order
            % of each node's position in the child list of its parent. Y
            % coordinates are (negative) depth in the tree.
            % Inputs:
            %   cvalfun
            %     Handle to a function that returns information about each
            %     node. Function must take a node object as input and
            %     return a scalar. Colors of the plotted nodes will
            %     represent the returned values. The data cursor can also
            %     be used to click on a node and read out its value from
            %     the Z coordinate.
            
            if nargin < 2
                cvalfun = [];
            end
            
            % info for each node
            ids = zeros(1, obj.nnodes);
                % id
            d = zeros(1, obj.nnodes);
                % depth
            x = zeros(1, obj.nnodes);
                % plotted x coordinate
            cvals = zeros(1, obj.nnodes);
                % color value
            
            % graph adjacency matrix
            A = false(obj.nnodes);
                % A(i,j) = true iff node with index j is a child of node
                % with index i
            
            % each node is given an index (into ids, d, x, cvals, A)
            % next index to use
            next_idx = 1;
            
            % x coordinate to use for the next leaf node we encounter
            % how much to increase it by for each leaf node
            next_leaf_x = 1;
            delta_x = 1;
            
            % recursively build ids, d, x, cvals, A
            get_info(obj.root);
            
            % plot figure
            gplot(A, [x(:), -d(:)]);
            set(get(gca, 'children'), 'color', [1 1 1] .* 0.5);
            hold_off = ~ishold;
            hold on;
            if isempty(cvalfun)
                h = scatter(x, -d, 'filled');
                set(h, 'MarkerFaceColor', [0, 0.447, 0.741]);
            else
                h = scatter(x, -d, [], cvals, 'filled');
                set(h, 'ZData', cvals);
                colorbar;
            end
            set(gca, 'ytick', -obj.depth : 0);
            grid on;
            xlim([min(x) - 1, max(x) + 1]);
            if hold_off
                hold off;
            end

            % recursively get node info
            function [ idx ] = get_info( node )
                % Inputs:
                %   node = current node (TreeNode)
                % Outputs:
                %   idx = index of of current node
                % Modifies vars in parent function:
                %   ids, d, x, cvals, A, next_idx, next_x, next_leaf_x
                
                % index to use for this node, next node
                idx = next_idx;
                next_idx = next_idx + 1;
                
                % set node id and depth
                ids(idx) = node.id;
                d(idx) = node.depth;
                
                % set color value if requested
                if ~isempty(cvalfun)
                    cvals(idx) = cvalfun(node);
                end
                
                % if this is a leaf node
                if isempty(node.children)
                    % set x coordinate, update it for next leaf node
                    x(idx) = next_leaf_x;
                    next_leaf_x = next_leaf_x + delta_x;
                
                else
                    % recursively call get_info() on children
                    % build a list of their indices
                    child_idx = zeros(1, numel(node.children));
                    for it = 1 : numel(node.children)
                        child_idx(it) = ...
                            get_info(node.children(it));
                    end
                    
                    % update adjacency matrix
                    A(idx, child_idx) = true;
                    
                    % center x coordinate over children
                    x(idx) = mean(x(child_idx));
                end
            end
        end % plot()

    end % methods
end % classdef