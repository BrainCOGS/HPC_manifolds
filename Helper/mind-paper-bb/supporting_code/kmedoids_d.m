% kmedoids_d()
% Use k-medoids to cluster dataa, given as a distance matrix
%
% Attempts to cluster data points such that the sum of distances between each
% point and the medoid of its assigned cluster is minimized. The medoid of a
% cluster is the data point with minimum average distance to all other points
% in its cluster.
%
% Usage:
%   [ as, m, cost ] = kmedoids_d( d, k, options, ... )
%
% Inputs:
%   d
%     A pairwise distance matrix. d(i, j) contains the distance between data
%     points i and j. d must be symmetric and contain zeros on the diagonal.
%   k
%     The number of clusters
%   options
%     A set of name/value pairs specifying optional parameters:
%     'init'
%       Explicit initialization conditions. Can be used to specify the initial
%       medoids as a length(k) vector m0, where m0(i) contains the index of the
%       data point serving as the medoid of cluster i. Entries of m0 must be
%       unique. Can also be used to specify the initial cluster asignments as a
%       length(n) vector as0, where as0(i) contains an integer from 1 to k
%       indicating the cluster assigned to data point i.
%       [Default: k random data points will be selected as the initial medoids]
%     'replicates'
%       Specifies the number of times to run the algorithm. On each run, k
%       random data points wil be selected as the initial medoids. The
%       clustering with the lowest cost (described below) will be returned.
%       Incompatible with the 'init' option. [Default: 1]
%     'verbose'
%       A boolean value specifying whether to print progress updates
%       [Default: false]
%
% Outputs:
%   as
%     A vector of cluster asignments. as(i) contains an integer from 1 to k
%     specifying the cluster assigned to data point i.
%   m
%     A vector of medoids. m(i) contains the index of the data point serving as
%     the medoid of cluster i.
%   cost
%     The sum of distances between each datapoint and the medoid of its
%     assigned cluster
%
% Notes:
%   Returns a local minimum solution. The 'replicates' option can be used to
%   search multiple local minima.

% Ver. 1.1.0, R. Low, 2015-06-02


function [ as, m, cost ] = kmedoids_d( d, k, varargin )

	% parse options into struct
	options = get_options(varargin);

	% primitive checks on check distance
	if numel(size(d)) ~= 2 || size(d, 1) ~= size(d, 2)
		error('Invalid distance matrix');
	end

	% initialize
	as = []; % cluster assignments
	m = []; % medoids
	cost = inf; % cost

	% loop thru replicates
	for rep = 1 : options.replicates

		% cluster data
		[my_as, my_m, my_cost] = do_kmedoids(d, k, options);

		% keep current clustering if its cost is the best so far
		if my_cost < cost
			as = my_as;
			m = my_m;
			cost = my_cost;
		end

		% print
		if options.verbose
			fprintf( ...
				'### iteration %d/%d: cost=%g\n', ...
				rep, options.replicates, my_cost ...
			);
		end
	end
end


function [ as, m, cost ] = do_kmedoids( d, k, options )
% Performs k-medoids clustering
%
% Inputs/outputs (as described in file header)
%  d: Pairwise distance matrix
%  k: Number of clusters
%  options: Options struct, as returned by get_options()
%  as: Vector of cluster asignments
%  m: Vector of medoids
%  cost: Cost of clustering
	
	% how many data points
	npts = size(d, 1);

	% initialize medoids
	% m(i) = index of data point that serves as the medoid of cluster i

	% default: select k random data points as initial medoids
	if isempty(options.init)
		idx = randperm(npts);
		m = idx(1 : k);

	% initialization specified
	else
		% initial medoids specified
		if numel(options.init) == k
			m = options.init;
		% initial cluster assignments specified; find corresponding medoids
		elseif numel(options.init) == npts
			%m = get_medoids(d, as, k);
			m = get_medoids(d, options.init, k);
		else
			error('Init option must specify medoids or cluster assignments');
		end
	end

	% loop until convergence
    k = numel(m);
	done = false;
	ctr = 1;
	while ~done
		if options.verbose
			fprintf('%d\n', ctr);
		end
		ctr = ctr + 1;

		% assign each data point to cluster with the closest medoid
 		[ as, cost ] = get_assignments( d, m );

		% find new medoids
		m_new = get_medoids( d, as, k );

		% terminate if medoids haven't changed
		if all(m_new == m)
			done = true;
		else
			m = m_new;
		end
	end
end


function [ m ] = get_medoids( d, as, k )
% Finds medoids, given clusters assigned to data points
% Inputs/outputs (as described in file header)
%  d: Pairwise distance matrix
%  as: Vector of cluster asignments
%  k: Number of clusters
%  m: Vector of medoids

	% loop thru clusters, find medoid of each
	m = zeros(1, k);
	for it = 1 : k

		% points in current cluster
		my_pts = find(as == it);
		if isempty(my_pts)
			%error('Empty cluster created');
            continue;
			% this shouldn't happen because a medoid is always a member of its
			% own cluster
		end

		% medoid is point with minimum distance to other points in cluster
		my_d = d(my_pts, my_pts);
		[ans, minloc] = min(sum(my_d, 2));
		m(it) = my_pts(minloc);
    end    
    m( m == 0 ) = [];    
end


function [ as, cost ] = get_assignments( d, m )
% Assigns data points to clusters given cluster medoids
% Inputs/outputs (as described in file header)
%  d: Pairwise distance matrix
%  m: Vector of medoids
%  as: Vector of cluster asignments
%  cost: Cost of returned clustering
		
		% assign each data point to cluster with the closest medoid
		[dist, as] = min(d(:, m), [], 2);
		% as(i) = index of cluster assigned to point i
		% dist(i) = distance from point i to closest medoid

		% sum of distances from each point to medoid of its assigned cluster
		cost = sum(dist);
end


function [ options ] = get_options( args )
% Parses options input to main function into a struct. Sets some default values.
% Inputs/outputs:
%   args: Cell array containing name/value pairs. Names are case insensitive.
%   options: Struct, where options.name = value. All field names lowercase.

	% is this a valid set of name/value pairs?
	if mod(numel(args), 2) ~= 0 || ~all(cellfun(@isstr, args(1 : 2 : end - 1)))
		error('Options must be specified as name/value pairs');
	end

	% loop thru name/value pairs, set corresponding struct fields/values
	options = struct;
	for it = 1 : 2 : numel(args) - 1
		try
			options.(lower(args{it})) = args{it + 1};
		catch
			error(sprintf('Error setting option ''%s''', args{it}));
		end
	end

	% cannot specify both init and replicates options
	if isfield(options, 'init') && isfield(options, 'replicates')
		error('The ''init'' and ''replicates'' options are incompatible');
	end

	% default: single replicate
	if ~isfield(options, 'replicates')
		options.replicates = 1;
	end

	% default: use random initialization
	if ~isfield(options, 'init')
		options.init = [];
	end

	% default: don't print output
	if ~isfield(options, 'verbose')
		options.verbose = false;
	end
end
