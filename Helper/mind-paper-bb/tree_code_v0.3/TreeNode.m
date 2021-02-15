% TreeNode class (abstract)
% Represents nodes of a decision tree.
%
% Inherits from handle class, so all copies reference the same underlying
% object.

classdef TreeNode < handle
    
    properties
        parent; % (TreeNode object)
            % Parent of this node. Empty if this is the root.
        children; % (Vector of TreeNode objects)
            % Children of this node. Empty if this is a leaf.
        id; % (integer)
            % A unique integer identifying this node.
            % Tree class will handle setting this property.
        depth; % (integer)
            % Depth of the node in the tree (root node has depth 0)
        train_pts; % (integer)
            % Number of points used to train node
        train_loss; % (scalar)
            % Value of loss function for local model on training points
        params; % (HandleVar object)
            % Parameters shared by all nodes in a tree should be stored as
            % fields of struct params.data. A HandleVar object is used
            % to avoid storing local copies for every node (HandleVar is a
            % handle class). Every new node in the tree must receive the
            % same params object.
            %
            % Must include fields:
            %   max_depth
            %     Maximum depth of the tree (root node has depth 0)
            %   min_leaf_pts;
            %     Minimum number of training points per leaf node.
            %   min_leaf_loss;
            %     A node can only be split if the loss of its local model
            %     exceeds this value. Corresponds to an average over local
            %     training points (e.g. entropy, MSE), not sum.
    end
    
    % Methods that must be implemented by subclasses
    methods (Abstract)
        
        % Create a new, empty node with same shared parameters
        spawn( obj );
            % Outputs:
            %   obj2
            %     An empty node. The shared 'params' HandleVar object
            %     is copied from the input to the new node.
            % Notes:
            %   The subclass must define this because there's currently no
            %   good way to call the subclass constructor here.
        
        % Train node's local model
        train( obj, x, y );
            % Inputs:
            %   x, y
            %     Matrix of input/output training data. Rows correspond to
            %     points, columns to dimensions.
            % Outputs:
            %   loss
            %     Value of the loss function for the local model on the
            %     training data. Corresponds to an average over points
            %     (e.g. entropy, MSE), not sum.
            % Postconditions:
            %   'train_pts', 'train_loss', and properties representing
            %   local model parameters are set.
        
        % Check whether a node can be split
        is_splittable( obj )
            % Outputs:
            %   splittable
            %     True if it's possible to split the given node. False if
            %     given node can't be split. This can happen if the node's
            %     loss is already below the minimum, the node contains too
            %     few training points to assign the minimum number to each
            %     child, or the maximum depth has already been hit.
        
        % Train node's gating function and create new child nodes with
        % trained local models
        split( obj, x, y );
            % Inputs:
            %   x, y
            %     Matrix of input/output training data. Rows correspond to
            %     points, columns to dimensions.
            % Outputs:
            %   as
            %     Assignment of training points to newly created child
            %     nodes. Logical matrix. as(i,j) is true if point i
            %     in training data is assigned to new child node j. Empty
            %     if split not successful (see notes).
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
            %   Split must respect the minimum number of points per node,
            %   specified by params.data.min_leaf_pts.
            %
            %   A split can fail if the new loss exceeds the loss for the
            %   current node. The new loss is given by the average loss
            %   over new child nodes, weighted by the number of points
            %   assigned to each.
        
        % Gating function. Assign data points to child nodesw
        assign( obj, x );
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions. Number of dimensions must match node's gating
            %     function.
            % Outputs:
            %   as
            %     Logical matrix specifying the child node to which each
            %     data point is assigned. as(i,j) = true if point i is
            %     assigned to child node j, otherwise false.
            % Preconditions:
            %   Gating function must be trained
        
        % Use local model to predict output, given input
        predict( obj, x );
            % Inputs:
            %   x
            %     Input data matrix. Rows correspond to points, columns to
            %     dimensions. Number of dimensions must match node's local
            %     model.
            % Outputs:
            %   yhat
            %     Matrix containing predicted output for each input. Rows
            %     correspond to points, columns to dimensions.
            % Preconditions:
            %   Node's local model must be trained.

    end % methods
end % classdef