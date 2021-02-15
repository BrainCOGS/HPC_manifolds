% greedyvq()
% Performs vector quantization using an iterative, greedy approach
%
% Constructs a Voronoi partitition of the input space using a set of code
% vectors, where every point is assigned to the nearest code vector (using
% l2 distance). Adds a new code vector at each step, chosen to be the data point
% with the greatest quantization error (distance to its currently assigned code
% vector. Termination occurs when the number of code vectors reaches a specified
% maximum or the goodness of fit exceeds a specified R^2 value.
%
% Usage:
%   [ as, c, dmin, r2 ] = greedyvq( x, kmax, options )
%
% Inputs:
%   x
%     Data points. Rows correspond to points, columns to dimensions
%   kmax
%     The maximum number of code vectors. Termination will occur when this
%     number is reached.
%
% Options (given as name/value pairs):
%   'c0'
%     An initial set of code vectors. Rows correspond to code vectors, columns
%     to dimensions. [Default: Pick a single data point uniformly at random]
%   'r2max'
%     Maximum goodenss of fit (R^2). Termination will occur if R^2 reaches
%     this value. [Default: 1.0]
%   'accelerate'
%     Boolean value. If true, accelerates the computation using the triangle
%     inequality to avoid unnecessary distance calculations. In cases where
%     distance calculations cannot be avoided (e.g. when the data are
%     intrinsically high dimensional and unstructured), this option adds a small
%     amount of extra overhead, and can be disabled. [Default: true]
%
% Outputs:
%   as
%     Code vector assignments. as(i) is an integer indicating the code vector 
%     assigned to data point i.
%   c
%     Code vectors. Rows correspond to code vectors, columns to dimensions.
%   dmin
%     Distance from each data point to its final, assigned code vector.
%   r2
%     A vector indicating the goodness of fit at each step of the procedure.
%     r2(k) contains the R^2 value when using code vectors 1 thru k.
%
% Notes:
%   Can use the triangle inequality to accelerate the computation by avoiding
%   unnecessary distance calculations, as described in [1] (lemma 1).
%
%   Memory usage: Let n be the number of points and d the number of dimensions.
%   This function stores O(kmax*(d + 1) + 2*n) elements of type double (not
%   including storage needed by external functions called). An extra kmax^2
%   elements are required when acceleration is enabled.
%
%   Because of the greedy, sequential nature of the algorithm, the code vectors
%   form a nested sequence. To examine multiple choices of k (for k < kmax)
%   take the first k code vectors.
%
%   Requires Mathworks statistics toolbox
%
% References:
%   [1] Elkan C. 2003. Using the triangle inequality to accelerate k-means


function [ as, c, dmin, r2 ] = greedyvq( x, kmax, varargin )

    % data size
    [npts, ndims] = size(x);

    % check kmax
    if isempty(kmax) || kmax == inf
        kmax = npts;
    end
    if isfinite(kmax) && (kmax < 1 || kmax > npts)
        error([ ...
            'Number of code vectors must be from 1 to the number of' ...
            ' data points' ...
        ]);
    end

    % parse options
    opt = get_options(x, kmax, varargin);

    % initialize

    % current number of code vectors
    k = size(opt.c0, 1);

    % code vectors
    c = zeros(kmax, ndims);
        % size(code vectors, dimensions)
    c(1 : k, :) = opt.c0;

    % code vector assignments
    [as, dmin] = knnsearch(opt.c0, x);
        % as(i) = index of code vector assigned to point i
        % dmin(i) = distance from point i to its assigned code vector

    % goodness of fit: r^2 for entire quantization
    varsum = sum(var(x, 1, 1));
    r2 = nan(kmax, 1);
    r2(k) = 1 - mean(dmin .^ 2) ./ varsum;
        % r2(i) contains r^2 for step i (using prototypes 1 thru i)

    % pairwise distances between code vectors
    % used to accelerate computation using triangle inequality, as in [1]
    if opt.accelerate
        % cd(i,j) = distance between code vectors i and j
        cd = zeros(kmax, kmax);
        if k > 1
            cd(1 : k, 1 : k) = squareform(pdist(opt.c0));
        end
    end

    % main loop

    % add prototypes until stopping criteria reached
    while k < kmax && r2(k) < opt.r2max
        k = k + 1;

        % select the point with the greatest error (furthest from currently
        % assigned code vector). add it as a new code vector.
        [maxval, maxloc] = max(dmin);
        c(k, :) = x(maxloc, :);

        % accelerate computation using the triangle inequality. only a subset
        % of points need be considered as candidates for reassigning to new
        % code vector. this avoids unnecessary distance computations.
        % use lemma 1 from [1]: if d(b,c) >= 2*d(x,b) then d(x,c) >= d(x,b)
        if opt.accelerate

            % update pairwise distances between code vectors
            cd(1 : k, k) = pdist2(c(1 : k, :), c(k, :));
            cd(k, 1 : k) = cd(1 : k, k)';

            % find points that are candidates for being reassigned to new code
            % vector
            candidates = find(cd(as, k) < 2 .* dmin);

            % a point is a candidate if the distance from its currently assigned
            % code vector to the new code vector is less than 2 times the
            % distance from the point to its currently assigned code vector

            % find distances from candidate points to new code vector
            dnew = pdist2(x(candidates, :), c(k, :));

            % select candidate points that are closer to new code vector than to
            % their previously assigned code vector
            sel = dnew < dmin(candidates);
                % sel is an indicator vector over candidates

            % reassign selected points to new code vector
            as(candidates(sel)) = k;
            dmin(candidates(sel)) = dnew(sel);
        
        % if no acceleration requested, all points are candidates
        % avoid book keeping overhead
        else

            % find distances from all points to new code vector
            dnew = pdist2(x, c(k, :));

            % select points that are closer to new code vector than to their
            % previously assigned code vector
            sel = dnew < dmin;

            % reassign selected points to new code vector
            as(sel) = k;
            dmin(sel) = dnew(sel);

        end % acceleration check

        % calculate new goodness of fit
        r2(k) = 1 - mean(dmin .^ 2) ./ varsum;

    end % main loop

    % truncate results to the number of code vectors actually used
    % may be less than kmax if r^2 stopping threshold reached
    c = c(1 : k, :);
    r2 = r2(1 : k);
end


function [ opt ] = get_options( x, kmax, args )
    % parses options, sets default values
    %
    % x, k = as passed to main function
    % args = cell array containing name/vlaue pairs, passed as varargin to
    %   main function
    % opt = options struct. field names are option names

    % create options struct
    if mod(numel(args), 2) ~= 0
        error('Options must be passed as name/value pairs');
    end
    try
        opt = struct(args{:});
    catch
        error('Invalid options');
    end

    % check options, set defaults

    % r squared
    if ~isfield(opt, 'r2max')
        opt.r2max = 1;
    end
    if opt.r2max == inf
        opt.r2max = 1;
    elseif opt.r2max > 1
        error('R^2 cannot exceed 1');
    end

    % triangle inequality acceleration
    if ~isfield(opt, 'accelerate')
        opt.accelerate = true;
    end

    % initial code vectors

    % if not specified, choose a single data point uniformly at random
    if ~isfield(opt, 'c0')
        opt.c0 = x(randi(size(x, 1)), :);

    % otherwise use specified values
    else
        if isempty(opt.c0)
            error('Initial code vectors are empty');
        end
        if size(opt.c0, 1) > kmax
            error('Number of initial code vectors exceeds kmax');
        end
        if size(opt.c0, 2) ~= size(x, 2)
            error('Initial code vectors must have same dimensionality as data');
        end
    end
end
