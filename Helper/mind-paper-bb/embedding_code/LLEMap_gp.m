classdef LLEMap_gp < handle
    
    properties
        Xtrain; % (matrix)
        % Input training points (points x dimensions)
        Ytrain; % (matrix)
        % Output training points (points x dimensions)
        Xdim; % (scalar)
        % Number of input dimensions
        Ydim; % (scalar)
        % Number of output dimensions
        model; % (cell array)
        % trained gp model (output dimensions x 1)
    end % properties
    
    methods
        
        % Constructor
        function [ obj ] = LLEMap_gp()
            %
        end % constructor
        
        
        % Fit mapping
        function [] = fit(obj, X, Y, options)
            % --------------------------------------------------------------------            
            % Inputs
            % --------------------------------------------------------------------
            %   X
            %     Matrix of training input points (points x dimensions)
            %   Y
            %     Matrix of training output points (points x dimensions)
            %   options
            %     Options for running fitrgp.
            
            if ~exist('options','var')
                options = [];
            end
            
            % --------------------------------------------------------------------
            % store training data
            % --------------------------------------------------------------------
            npts = size(X, 1);
            if size(Y, 1) ~= npts
                error([ ...
                    'Inputs and outputs must contain same number of ' ...
                    'points' ...
                    ]);
            end
            obj.Xtrain = X;
            obj.Ytrain = Y;
            Xdim = size(X, 2);
            Ydim = size(Y, 2);
            obj.Xdim = Xdim;
            obj.Ydim = Ydim;
            
            % --------------------------------------------------------------------
            % get parameter options for training
            % --------------------------------------------------------------------
            if isfield(options, 'FitMethod')
                FitMethod = options.FitMethod;
            else
                if size(X, 1) <= 2000
                    FitMethod = 'exact';
                else
                    FitMethod = 'sd';
                end
            end
            if isfield(options, 'ActiveSetSize')
                ActiveSetSize = options.ActiveSetSize;
                if ActiveSetSize <= 1
                    ActiveSetSize = round(ActiveSetSize * size(X,1));
                end
                ActiveSetSize = min(size(X, 1), ActiveSetSize);
            else
                ActiveSetSize = min(size(X, 1), 2000);
            end
            if isfield(options, 'ActiveSetMethod')
                ActiveSetMethod = options.ActiveSetMethod;
            else
                ActiveSetMethod = 'random';
            end
            if isfield(options, 'NumActiveSetRepeats')
                NumActiveSetRepeats = options.NumActiveSetRepeats;
            else
               NumActiveSetRepeats = 3;
            end            
            
            % --------------------------------------------------------------------
            % train model
            % --------------------------------------------------------------------
            model = cell(Ydim,1);
            parfor id = 1:Ydim
                fprintf(...
                    'Training gpr model on output dimension %d of %d\n', ...
                    id, Ydim);
                thisY = Y(:,id);                
                model{id} = fitrgp(X, thisY, 'FitMethod', FitMethod, ...
                    'ActiveSetSize', ActiveSetSize, 'ActiveSetMethod', ...
                    ActiveSetMethod, 'NumActiveSetRepeats', NumActiveSetRepeats);
                model{id} = compact(model{id});
            end
            obj.model = model;
        end % fit()
        
        % Map new input points into output space
        function [ Y ] = transform(obj, X)
            
            % check inputs
            if isempty(obj.Xtrain)
                error('Model must be fit first');
            end
            if size(X, 2) ~= obj.Xdim
                error('Dimensionality of points doesn''t match model');
            end
            
            Ydim = obj.Ydim;
            Y = zeros(size(X, 1), Ydim);
            for id = 1:Ydim
                Y(:, id) = obj.model{id}.predict(X);
            end
        end % transform()
        
    end % methods
end % classdef