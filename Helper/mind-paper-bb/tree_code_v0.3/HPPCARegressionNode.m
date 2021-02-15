% HPPCARegressionNode class
% Hybrid PPCA regression. Same as PPCARegressionNode, but uses a hybrid
% training procedure for speed when optimizing splits. For each candidate
% splitting direction, optimizes the splitting threshold assuming spherical
% Gaussian distributions. Once the threshold has been determined, fits PPCA
% models to data on each side of the split. The value of the loss function
% for the PPCA models is then used to pick the best direction.

classdef HPPCARegressionNode < PPCARegressionNode
    
    methods
        
        % Constructor
        function [ obj ] = HPPCARegressionNode( varargin )
            % just call the superclass constructor
            obj = obj@PPCARegressionNode(varargin{:});
        end
        
        
        % Create a new, empty node with same shared parameters
        function [ obj2 ] = spawn( obj )
            % Outputs:
            %   obj2
            %     An empty node. Its 'params' HandleVar object is copied
            %     from obj.
            
            obj2 = HPPCARegressionNode([]);
            obj2.params = obj.params;
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
            %   Chooses the splitting threshold to minimize the loss,
            %   assuming spherical Gaussian distributions. Searches over
            %   all splitting threshold indices using fast procedure based
            %   on running sums. The spherical assumption is dropped after
            %   finding the splitting thresholds. PPCA models are fit to
            %   the data in each candidate child node, and used to compute
            %   the loss for the split.
            
            % number of points, output dimensions
            [npts, nd] = size(y);
            
            % min/max threshold index such no child node receives fewer
            % than the minimum number of points (this could be violated if
            % many points share the same projections)
            t_min = obj.params.data.min_leaf_pts;
            t_max = npts - t_min;
                % t is an index of the sorted projections s.t. all points
                % w/ index > t are assigned to the right node
                
            % find threshold index using fast procedure from RegressionNode
            
            % how many points in left/right node, for each threshold
            nL = (t_min : t_max)';
            nR = npts - nL;
            
            % sum of output vectors in left/right node
            % for each threshold
            csY = cumsum(y, 1);
            sL = csY(t_min : t_max, :);
            sR = bsxfun(@minus, csY(end, :), sL);
            
            % value of surrogate loss function for each threshold
            surloss = -sum(sL.^2, 2) ./ nL - sum(sR.^2, 2) ./ nR;
                % minimizing this loss function is equivalent to minimizing
                % the average mse in left/right nodes, weighted  by the
                % number of points assigned to each. minimizing the squared
                % error is equivalent to maximizing the likelihood of
                % spherican gaussians.
            
            % find best threshold index
            [~, minloc] = min(surloss);
            t = nL(minloc);
            
            % fit ppca models
            loss_L = nodes(1).train([], y(1:t, :));
            loss_R = nodes(2).train([], y(t+1:end, :));
            child_loss = [loss_L, loss_R];
            
            % loss for the split, based on ppca models
            npts_L = t;
            npts_R = npts - t;
            loss = npts_L/npts*loss_L + npts_R/npts*loss_R;
        end % split_opt_1d()
    
    end % methods
end % classdef