% BSHNode class (abstract)
% Implements nodes of a decision tree using hyperplanes to partition the
% input space.

% Performs binary splits using axis-parallel or arbitrarily-oriented
% hyperplanes. Can perform random splits, or search over a randomly chosen
% set of splitting directions. In the latter case, the splitting threshold
% along each direction is optimized, and the split that most reduces the
% loss function is chosen.


classdef BSHNode < TreeNode
    
    properties
        
        % Split parameters
        split_dir; % (scalar or column vector)
            % If using axis parallel splits, split_dir is a scalar
            % specifying the dimension to split along. Otherwise, split_dir
            % is the normal vector of the splitting hyperplane.
        split_thresh; % (scalar)
            % Splitting threshold. Data points whose projection onto the
            % split direction is <= this value will be assigned to child
            % node 1, otherwise child node 2.
            
        % Shared parameters stored in params object:
        % (must be set by subclass constructor)
        %
        % ndir
        %   Number of split directions to try when optimizing splits.
        % axis_parallel
        %   If true, splitting hyperplanes will be axis parallel. If false,
        %   they'll be randomly oriented.
        % random_split
        %   If true, a random threshold and direction will be chosen for
        %   each split. Otherwise, they'll be optimized.
        %
        % The following shared parameters are also required
        % (from TreeNode):
        %
        % max_depth
        %   Maximum depth of the tree (root node has depth 0)
        % min_leaf_pts
        %   Minimum number of training points per leaf node
        % min_leaf_loss
        %   A node can only be split if the loss of its local model exceeds
        %   this value. Corresponds to an average over local training
        %   points (e.g. entropy, MSE), not sum.

    end % properties
    
    
    methods (Abstract)
        
        % Subclass must also implement abstract TreeNode methods that
        % aren't implemented here: train(), predict()
        
        % Find optimal splitting threshold along a single direction
        split_opt_1d( obj, x, y, nodes );
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
            %   Chooses the splitting threshold to minimize the loss. Split
            %   must respect the minimum number of points per node
            %   specified by params.data.min_leaf_pts. The loss need not
            %   obey any constraints.
            %
            %   Node objects are passed as input for speed. This allows
            %   recycling existing objects instead of constructing new
            %   ones. Properties may contain existing values from previous
            %   usage.
            
    end % abstract methods
    
    
    methods
        
        % Check whether a node can be split
        function [ splittable ] = is_splittable( obj )
            % Outputs:
            %   splittable
            %     True if it's possible to split the given node, otherwise
            %     false.
            % Notes:
            %   A given node can be split iff its depth is below the
            %   maximum, its loss is above the minimum, and it contains
            %   enough training points such that new child nodes receive at
            %   least the minimum number.
            
            splittable = ( ...
                obj.depth < obj.params.data.max_depth ...
                && obj.train_pts >= 2 * obj.params.data.min_leaf_pts ...
                && obj.train_loss > obj.params.data.min_leaf_loss ...
            );
        end % is_splittable()
            
        
        % Split a node to produce new child nodes
        function [ as ] = split( obj, x, y )
            % Inputs:
            %   x, y
            %     Input/output training points assigned to given node. Rows
            %     correspond to points, columns to dimensions.
            % Outputs:
            %   as
            %     Assignment of training points to newly created child
            %     nodes. Logical matrix. as(i,j) = true iff point i in
            %     training data is assigned to new child node j. Empty
            %     if split not successful.
            % Preconditions:
            %   The node contains a trained local model and no children.
            %   The node is splittable, as determined by is_splittable().
            % Postconditions:
            %   If the split is successful, the current node's gating
            %   function is trained, and new child nodes are created. Each
            %   child node contains a trained local model, and its 'depth',
            %   'train_pts', 'train_loss', and 'params' properties are set.
            %   Child nodes are linked into the tree by setting their
            %   'parent' property to the current node, and the current
            %   node's 'children' property to a vector of children. If the
            %   split is unsuccessful, no changes are made.
            % Notes:
            %   A split is successful if the average loss over new child
            %   nodes (weighted by the number of points in each) is less
            %   than the loss for the current node.
            
            % check that node has no children
            if ~isempty(obj.children)
                error('Can''t split a non-leaf node');
            end
            
            % training data size
            npts = size(x, 1);
            if size(y, 1) ~= npts
                error('x and y must contain equal numbers of points');
            end
            
            % find the best split
            if obj.params.data.random_split
                [ split_dir, split_thresh, idx, nodes, loss ] = ...
                    obj.split_rand( x, y );
            else
                [ split_dir, split_thresh, idx, nodes, loss ] = ...
                    obj.split_opt( x, y );
            end
            
            % perform split if we reduced the loss sufficiently
            if loss < obj.train_loss
                success = true;
                
                % save split parameters
                obj.split_dir = split_dir;
                obj.split_thresh = split_thresh;
                
                % link child nodes into tree
                nodes(1).parent = obj;
                nodes(2).parent = obj;
                child_depth = obj.depth + 1;
                nodes(1).depth = child_depth;
                nodes(2).depth = child_depth;
                obj.children = nodes;
                
                % return assignment of training points to children
                as = [~idx, idx];
            
            % if split not successful, return null result
            else
                as = [];
            end
        end % split()
        
        
        % Find a good split by searching over multiple splitting
        % directions and optimizing the splitting threshold for each
        function [ split_dir, split_thresh, idx, nodes, loss ] = ...
                split_opt( obj, x, y )
            % Inputs:
            %   x, y
            %     Input, output points. Rows correspond to points, columns
            %     to dimensions
            % Outputs:
            %   split_dir
            %     Split direction. If using axis parallel splits, given as
            %     an integer specifying the axis to split along. Otherwise,
            %     given as the normal vector of the splitting hyperplane.
            %   split_thresh
            %     Splitting threshold. Points are assigned to the right
            %     child node if their projection onto the split direction
            %     exceeds threshold, otherwise to the left.
            %   idx
            %     Binary assignment vector. idx(i) = true if point i is
            %     asigned to the right child node, otherwise 0.
            %   nodes
            %     Left/right child nodes produced by the split. Given as
            %     a vector [node_L, node_R] of node objects with trained
            %     local models. Set properties also include 'train_pts',
            %     and 'train_loss'. Child nodes aren't linked into the tree
            %     (the 'parent' property isn't set).
            %   loss
            %     Value of loss function for the computed split. Given by
            %     the average loss for the local model in the left/right
            %     child nodes (weighted by the number of points in each).
            % Notes:
            %   Chooses a random set of candidate split directions. Chooses
            %   the splitting threshold along each direction to minimize
            %   the loss function. Returns the split with minimum loss.

            % number of points, input dimensions
            [npts, nd] = size(x);
            
            % choose split directions, project data onto them
            
            % axis-parallel splits
            if obj.params.data.axis_parallel
                ndir = min(obj.params.data.ndir, nd);
                all_split_dir = randperm(nd, ndir);
                all_proj = x(:, all_split_dir);
            % randomly oriented splits
            else
                ndir = obj.params.data.ndir;
                all_split_dir = randn(nd, ndir);
                all_split_dir = bsxfun( ...
                    @rdivide, ...
                    all_split_dir, ...
                    sqrt(sum(all_split_dir .^ 2, 1)) ...
                ); % normalize
                all_proj = x * all_split_dir;
            end
            % all_split_dir(:, i) is the ith splitting direction
            % all_proj(:, i) containts the projections of all input points
            % onto it.
            
            % sort all projections in ascending order
            [all_proj, all_ord] = sort(all_proj, 1);
            
            % best split found so far
            loss = inf;
            nodes = [obj.spawn(), obj.spawn()];
            child_loss = [];
            di = 0;
                % index of best splitting direction
            t = 0;
                % threshold index. points 1:t assigned to left node,
                % points t+1:npts assigned to right node.
            
            % always work with pre-created nodes to avoid repeatedly
            % constructing new node objects
            my_nodes = [obj.spawn(), obj.spawn()];
            
            % search split directions, find optimal threshold for each
            for it = 1 : ndir
                
                % sort data in order of input projections along current
                % direction
                my_proj = all_proj(:, it);
                ord = all_ord(:, it);
                my_x = x(ord, :);
                my_y = y(ord, :);
                
                % get best split along current direction
                
                [my_t, my_child_loss, my_loss] = ...
                    obj.split_opt_1d( my_x, my_y, my_nodes );
                
                % keep split if it lowers the loss
                if my_loss < loss
                    di = it;
                    t = my_t;
                    loss = my_loss;
                    child_loss = my_child_loss;
                    
                    % swap best nodes and working nodes to reuse existing
                    % node objects and avoid constructing new ones
                    tmp = nodes;
                    nodes = my_nodes;
                    my_nodes = tmp;
                end
            end % loop over split directions
            
            % best split direction
            split_dir = all_split_dir(:, di);
            
            % choose threshold randomly between projection corresponding to
            % t and the next largest projection
            split_thresh = all_proj(t, di) + ...
                rand * (all_proj(t+1, di) - all_proj(t, di));
            
            % assign points to right node if projection > threshold
            idx = false(npts, 1);
            idx(all_ord(t+1 : end, di)) = true;
            
            % set number of training points and loss for each node
            nodes(1).train_pts = t;
            nodes(2).train_pts = npts - t;
            nodes(1).train_loss = child_loss(1);
            nodes(2).train_loss = child_loss(2);
        end % split_opt()
        
        
        % Choose random split direction and threshold
        function [ split_dir, split_thresh, idx, nodes, loss ] = ...
                split_rand( obj, x, y )
            % Inputs:
            %   x, y
            %     Input, output points. Rows correspond to points, columns
            %     to dimensions
            % Outputs:
            %   split_dir
            %     Split direction. If using axis parallel splits, given as
            %     an integer specifying the axis to split along. Otherwise,
            %     given as the normal vector of the splitting hyperplane.
            %   split_thresh
            %     Splitting threshold. Points are assigned to the right
            %     child node if their projection onto the normal vector
            %     exceeds threshold, otherwise to the left.
            %   idx
            %     Binary assignment vector. idx(i) = true if point i is
            %     asigned to the right child node, otherwise 0.
            %   nodes
            %     Left/right child nodes produced by the split. Given as
            %     a vector [node_L, node_R] of node objects with trained
            %     local models. Set properties also include 'train_pts',
            %     and 'train_loss'. Child nodes aren't linked into the tree
            %     (the 'parent' property isn't set).
            %   loss
            %     Value of loss function for the computed split. Given by
            %     the average loss for the local model in the left/right
            %     child nodes, weighted by the number of training points
            %     assigned to each.
            
            % number of points, input dimensions
            [npts, nd] = size(x);
            
            % choose split direction, project data onto it
            
            % axis-parallel split
            if obj.params.data.axis_parallel
                split_dir = randi(nd);
                proj = x(:, split_dir);
            % randomly oriented split
            else
                split_dir = randn(nd, 1);
                split_dir = split_dir ./ norm(split_dir);
                proj = x * split_dir;
            end
            % split_dir is the splitting direction
            % proj contains projections of the data onto it
            
            % choose splitting threshold from uniform distribution on an
            % interval such that no child node receives fewer than
            % minimum number of points (this could be violated if many
            % points share the same projections)
            sproj = sort(proj);
            a = obj.params.data.min_leaf_pts;
            b = npts - a + 1;
            split_thresh = sproj(a) + rand .* (sproj(b) - sproj(a));
            
            % assign training points to child nodes
            idx = proj > split_thresh;
            
            % create left/right child nodes
            % w/ copy of parent's params object
            node_L = obj.spawn();
            node_R = obj.spawn();
            nodes = [node_L, node_R];
            
            % train left/right child nodes
            loss_L = node_L.train(x(~idx, :), y(~idx, :));
            loss_R = node_R.train(x(idx, :), y(idx, :));
            
            % average loss, weighted by number of points in each child
            npts_R = nnz(idx);
            npts_L = npts - npts_R;
            loss = npts_R/npts * loss_R + npts_L/npts * loss_L;
        end % split_rand()
        
        
        % Assign data points to child nodes
        function [ as ] = assign( obj, x )
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions. Number of dimensions must match node's gating
            %     function.
            % Outputs:
            %   as
            %     Sparse logical matrix specifying the child node to which
            %     each data point is assigned. as(i,j) = true if point i
            %     is assigned to child node j, otherwise false.
            % Preconditions:
            %   Split parameters must be trained (node can't be a leaf
            %   node).
            % Notes:
            %   A point is assigned to child node 2 if its projection onto
            %   the normal vector of the splitting hyperplane exceeds the
            %   splitting threshold, otherwise child node 1.
            
            if obj.params.data.axis_parallel
                idx = x(:, obj.split_dir) > obj.split_thresh;
            else
                idx = x * obj.split_dir > obj.split_thresh;
            end
            as = [~idx, idx];
        end % assign()
        
    end % methods
end % classdef