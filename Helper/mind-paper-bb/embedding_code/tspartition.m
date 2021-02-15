% tspartition class
% Partition a time series into blocks for k fold cross validation.
%
% The time series is partitioned into k continuous blocks of equal length.
% Each block is used as a test set. The corresponding training set consists
% of the remaining time points. Training set points within a specified
% margin of the test set may optionally be omitted. This feature can be
% used to ensure that there are no temporal correlations between the
% training and test sets.
%
% Other than the constructor, the methods and properties of this class are
% identical to the Mathworks cvpartition class.

classdef tspartition
    
    % Hidden properties
    properties (GetAccess = public, SetAccess = protected, Hidden=true)
        
        % training/test set blocks
        test_start;
        test_end;
            % start/end index of each test set
        train_l_start;
        train_l_end;
            % start/end index of left part of each training set
        train_r_start;
        train_r_end;
            % start/end index of right part of each training set
        
        % defined by cvpartition:
        Type = 'kfold';
            % type of cross validation
        N;
            % how many data points
    end % hidden properties
    
    % Visible properties
    properties (GetAccess = public, SetAccess = protected)
        
        % defined by cvpartition:
        NumObservations;
            % how many data points
        NumTestSets;
            % how many test sets
        TrainSize;
            % how many training sets
        TestSize;
            % how many test sets
    end % visible properties
    
    
    methods
        
        % Constructor
        function [ obj ] = tspartition( t, k, margin )
            % Inputs:
            %   t
            %     A vector containing the time of each data point. Must be
            %     sorted in ascending order.
            %   k
            %     Number of folds for k-fold cross validation
            %   margin (optional)
            %     The minimum interval (in units of time) between the edges
            %     of a test set and its corresponding training set.
            
            t = t(:)';
            n = numel(t);
            
            % call parent constructor
            %obj = obj@cvpartition(n, 'kfold', k);
            
            % no margin by default
            if nargin < 3 || isempty(margin)
                margin = 0;
            end
            
            % divide time interval into k equal blocks (test sets)
            edges = interp1( ...
                t, 1:n, ...
                linspace(t(1), t(end), k+1), ...
                'nearest' ...
            );
            edges(end) = n+1;
            obj.test_start = edges(1 : end-1);
            obj.test_end = edges(2 : end) - 1;
            
            % start/end indices for left side of training set
            obj.train_l_start = ones(1, k);
            if margin == 0
                obj.train_l_end = obj.test_start - 1;
            else
                obj.train_l_end = zeros(1, k);
                for it = 2 : k
                    idx = find( ...
                        t < t(obj.test_start(it)) - margin, ...
                        1, 'last' ...
                    );
                    if ~isempty(idx)
                        obj.train_l_end(it) = idx;
                    end
                end
            end
            
            % start/end indices for right side of training set
            if margin == 0
                obj.train_r_start = obj.test_end + 1;
            else
                obj.train_r_start = zeros(1, k) + n + 1;
                for it = 1 : k-1                    
                    idx = find(t > t(obj.test_end(it)) + margin, 1);
                    if ~isempty(idx)
                        obj.train_r_start(it) = idx;
                    end
                end
            end
            obj.train_r_end = zeros(1, k) + n;
            
            obj.NumTestSets = k;
            obj.TrainSize = ( ...
                max(0, obj.train_l_end - obj.train_l_start) ...
                + max(0, obj.train_r_end - obj.train_r_start) ...
            ) + 2;
            obj.TestSize = obj.test_end - obj.test_start + 1;
            obj.NumObservations = n;
            obj.N = n;
        end % constructor
        
        
        % Rerandomize cross-validation partition
        % Only defined for compatibility with cvpartition
        function [] = repartition( obj )
            error('Can''t repartition a time series');
        end % repartition()
        
        
        % Get training set
        function trainIndices = training( obj, c )
            % Inputs:
            %   c
            %     Index (1-k) of the requested training set
            % Outputs:
            %   trainIndices
            %     A logical vector. trainIndices(i) = true if point i is
            %     in the cth training set, otherwise false.
            
            trainIndices = false(obj.N, 1);
            trainIndices(obj.train_l_start(c) : obj.train_l_end(c)) = true;
            trainIndices(obj.train_r_start(c) : obj.train_r_end(c)) = true;
        end % training()
        
        
        % Get test set
        function testIndices = test( obj, c )
            % Inputs:
            %   c
            %     Index (1-k) of the requested test set
            % Outputs:
            %   testIndices
            %     A logical vector. testIndices(i) = true if point i is
            %     in the cth test set, otherwise false.
            
            testIndices = false(obj.N, 1);
            testIndices(obj.test_start(c) : obj.test_end(c)) = true;
        end % test()
    
    end % methods
end % classdef