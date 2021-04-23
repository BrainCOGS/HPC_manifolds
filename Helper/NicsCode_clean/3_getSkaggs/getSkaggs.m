function out = getSkaggs(...
    input1DFFthresholding, input2log, input3gaussianSmoothingSigma, ...
    input4trialsAverageMap, input5numberShuffles, input6shuffleFilename, ...
    input7calcRealData, input8imagingFilename, input9thresholdSig, ...
    input10dimensions, input11binEdges, input12randomDimMethod, input13randomBinEdges, ...
    input14taskType, input15whichROIs, input16gaussianSmoothingToggle)
%
% calculates the mutual information and neuronal receptive fields along
% specified behavioral dimensions.
%
% example:
% out = getSkaggs([0 0], 'noLog', 0, 'keepTrials', 10, 'E47_ExY_%d.mat', 1, 'E47_20170927_70per_userSetSD5minDur0.modeling.mat', 2, {'Evidence', 'Position'}, {[-20:20], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
%
% *** INPUTS ***
%
% input1DFFthresholding - array of two integers, first integer sets the
% windowsize and second integer sets the thresholding
% if [0 0] == nothing happens
% if [0 4] == no filtering, threshold data 4*robustSTD
% if [11 4]== filter with windowsize=11 and threshold data 4*robustSTD
% if [11 0]== filter with windowsize=11, but don't threshold
%
% input2log - character vector
% 'log10' - take log10 of activity values with log10(v+1);
% 'noLog' - do not take log10 activity
%
% input3gaussianSmoothingSigma - double
% if ==0, no smoothing
% if >0, use value as sigma for gaussian smoothing
%
% input4trialsAverageMap - vector of integers, or character vector 'all'
% specify which trials should be used to construct the average maps for
% each neuron; useful for comparing correct vs error trials or left vs
% right trials. if specified as 'keepTrials', then construct average maps
% from only mainTrial and goodQuality trials from score.trial
%
%   ex:
%   'keepTrials' - do trials that are true for isGood and mainMaze
%   'keepTrials + leftChoice' - do 'keepTrials' and left choice trials
%   'keepTrials + rightChoice' - do 'keepTrials' and right choice trials
%   'goodTrials' - do trials that are true for isGood
%
% input5numberShuffles - integer 0 or >0
% specify number of shuffles to retrieve/write, using circshift within
% trial
% if == 0, do not do shuffle tests
% if >0, do that number of shuffles
%
% input6shuffleFilename - character vector
% character vector to either load or write shuffle tests to
% passed to sprintf(); needed in the following form with the '%d.mat' where
% the shuffle test number replaces %d
% Example: 'getSkaggs_shuffle%d.mat'
%
% input7calcRealData - integer 0 or 1
% if 0, do not calculate real data mutual information metric or average maps; write new
% shuffles
% if 1, calculate real data mutual information metric and average maps
%
% input8imagingFilename - character vector filename for animal session
% example: 'E22_20170215_30p.modeling_NEW.mat'
%
% input9thresholdSig - threshold for significance of comparing mutual
% information values with shuffle test values, e.g. 2 for 2SD above mean - ALSO
% used for the "shroud", i.e. threshold map, in pixelwise
%
% input10dimensions - cell array of character vectors specifying which
% dimension(s) to test for mutual information metric. case sensitive and must be found
% in output table of extractVariables()
% example: {'Evidence', 'Position', 'Choice'}
%
% input11binEdges - cell array of bin edges passed to discretize().
% ex: {[-20:20], [0:10:300]} for an ExY analysis
%
% input12randomDimMethod - 1-by-2 cell array of character vectors to
% specify which dimension to randomize with respect to a second dimension.
% Dim{1} is randomized with respect to its joint distribution with Dim{2}.
% ex: {'Evidence', 'Position'}  ...which is interpreted as 'randomly sample
% Evidence with respect to Position'
%
% input13randomBinEdges - 1-by-2 cell array of double vector. This sets the
% bin edges for the random dimensions in the same fashion of how
% input11binEdges sets the binEdges for the real dimension(s)
%
% input14taskType - 'towers' or 'alternation', used for extractVariables().
%
% input15whichROIs - which ROIs to use for mutual information metric and average map
% calculations. Specify either 'all' or a double vector of ROI labels
%   'all' - use all ROIs in modeling.mat file
%   double vector - use only ROIs in vector
%
% input16gaussianSmoothingToggle - character vector to specify which
% part(s) of the analysis to use gaussian filtering on:
%       - 'skaggsOnly' allows smoothing on skaggs calculation and sequence
%       plot
%       - 'pixelwiseOnly' allows smoothing on bin-wise analysis of mean DFF
%       in pixelwise
%       - 'both' allows smoothing on both analyses above
%
%
% *** OUTPUTS ***
%
% out.unfilteredAverageMaps - contains mean DFF (lambda(x)) for ROIs
% out.pixelwise - bin-wise comparison of real mean DFF to shuffled mean DFF
% out.pX - number of frames the animal spent in each bin
% out.skaggsMetric - summarizes mutual information analysis
% out.sequencePlot - creates a sequence plot with informative ROIs
% out.config - configuration settings that tend not to change
% out.argins - user-specified input arguments
% out.DFF - DFF data used to construct (lambda(x))
% out.dimensionIndices - indices used to index DFF data into bins
% out.trials - trial IDs used in the analysis
% out.ROIs - ROI IDs used in the analysis
%
%% Set configurations
tic;

config.circshiftIndLimit            = 1;    % limit for within-trial shuffling
config.shuffleRNG                   = 2018; % RNG set to this each time
if strcmp(input14taskType, 'towers')
    config.mazeRegion               = 2;    % cue + delay region
elseif strcmp(input14taskType, 'alternationJeff')
    config.mazeRegion               = 6;    
end
config.activityType                 = 2;    % DFF activity
config.trialLength                  = 0;    % sample DFF by frame
config.eventThresholdCoefficient    = 5;    % to detect events in DFF
config.discretizeDefaultNumberBins  = 10;   % default if not specified
config.gaussianSmoothingParameters  = 5;    % hsize in fspecial('gaussian')
config.shroudThreshold              = 2;    % number of SD to build the shroud at
config.random_pX_caxisLim           = 50;   % colorbar limit for randomized pX
config.seqPlotPvalueThreshold       = 0.05; % alpha level for sequence plot threshold
config.plotRandomizedVariable       = 0;    % if 1, make plot, otherwise don't
%% Save input arguments
argins.DFFthresholding              = input1DFFthresholding;
argins.log                          = input2log;
argins.gaussianSmoothingSigma       = input3gaussianSmoothingSigma;
argins.trialsAverageMap             = input4trialsAverageMap;
argins.numberShuffles               = input5numberShuffles;
argins.shuffleFilename              = input6shuffleFilename;
argins.calcRealData                 = input7calcRealData;
argins.imagingFilename              = input8imagingFilename;
argins.significanceThreshold        = input9thresholdSig;
argins.dimensions                   = input10dimensions;
argins.randomDimMethod              = input12randomDimMethod;
argins.taskType                     = input14taskType;
argins.whichROIs                    = input15whichROIs;
argins.binEdges                     = input11binEdges;
argins.randomBinEdges               = input13randomBinEdges;
argins.gaussianSmoothingToggle      = input16gaussianSmoothingToggle;

clear input*
%% Check inputs
helper_getSkaggs_checkInputs(argins);
%% Preallocate output variables
unfilteredAverageMaps = []; % contains lambda(x) for each ROI
pixelwise = []; % contains smoothed lambda(x) for real and shuffle data
pX = []; % the number of frames that the animal spends in each bin
skaggsMetric = []; % contains MI results
sequencePlot = []; % contains 1D sequence plot
%% Prepare data
wb = waitbar(0, 'Preparing data.');
% get data from extractVariables
tbl = extractVariables(argins.whichROIs, config.mazeRegion, ...
    argins.trialsAverageMap, config.activityType, config.trialLength, ...
    0, config.eventThresholdCoefficient, argins.DFFthresholding, ...
    argins.imagingFilename, 'none', argins.taskType, [], []);

data_behav = tbl.behavioralVariables;       % behavioral data
data_ROIs  = tbl.ROIactivities;             % ROI DFF data

numDimensions = length(argins.dimensions);  % number of dimensions
numROIs = size(data_ROIs, 2);               % number of ROIs

instTrialsAverageMap = unique(data_behav.Trial); % trial IDs
out.trials = instTrialsAverageMap;          % save trial IDs

clear D tbl
%% Extract the dimensions indices

dimInds = zeros(size(data_ROIs, 1), numDimensions); % keeps dimension indices
mapSize = []; % number of bins in each dimension
for i = 1:numDimensions % loop through each dimension
    % if the dimension is 'Random', then look at argins.randomDimMethod to
    % see which dimension to randomize. if the dimension is not 'Random',
    % then discretize the data from data_behav.
    if strcmp(argins.dimensions{i}, 'Random') == 0 && strcmp(argins.dimensions{i}, 'Random3D') == 0
        % if non-randomized dimension, get values from data_behav and
        % discretize.
        vals = table2array(data_behav(:, argins.dimensions{i}));
        dimInds(:,i) = helper_getSkaggs_discretize(vals, argins, config, i, 0);
        
        % calculate the maximum map size. for N bin edges, there are (N-1) bins
        mapSize(i) = length(argins.binEdges{i}) - 1;
    elseif strcmp(argins.dimensions{i}, 'Random') == 1
        % discretize dimension 1
        vals_dim1 = table2array(data_behav(:, argins.randomDimMethod{1}));
        vals_dim1 = helper_getSkaggs_discretize(vals_dim1, argins, config, i, 1);
        
        % discretize dimension 2
        vals_dim2 = table2array(data_behav(:, argins.randomDimMethod{2}));
        vals_dim2 = helper_getSkaggs_discretize(vals_dim2, argins, config, i, 2);
        
        % sample from dimension 1 with respect to dimension 2
        dimInds(:, i) = helper_getSkaggs_randomDimension(vals_dim1, vals_dim2);
        
        % plot the pX comparing the non-randomized and randomized dimension
        if config.plotRandomizedVariable==1
            plot_randomVariables_pX(vals_dim1, vals_dim2, dimInds(:,i), argins, config);
        end
        
        % calculate the maximum map size. for the randomized dimension, the
        % number of bins for this dimension is found using the first set of
        % argins.randomBinEdges.
        mapSize(i) = length(argins.randomBinEdges{1}) - 1;
    elseif strcmp(argins.dimensions{i}, 'Random3D') == 1
        % discretize dimension 1
        vals_dim1 = table2array(data_behav(:, argins.randomDimMethod{1}));
        vals_dim1 = helper_getSkaggs_discretize(vals_dim1, argins, config, i, 1);
        
        % discretize dimension 2
        vals_dim2 = table2array(data_behav(:, argins.randomDimMethod{2}));
        vals_dim2 = helper_getSkaggs_discretize(vals_dim2, argins, config, i, 2);
        
        % discretize dimension 3
        vals_dim3 = table2array(data_behav(:, argins.randomDimMethod{3}));
        vals_dim3 = helper_getSkaggs_discretize(vals_dim3, argins, config, i, 3);
        
        % sample from dimension 1 with respect to dimension 2
        dimInds(:, i) = helper_getSkaggs_randomDimension3D(vals_dim1, vals_dim2, vals_dim3);
        % plot the pX comparing the non-randomized and randomized dimension
        % plot_randomVariables_pX(vals_dim1, vals_dim2, dimInds(:,i), argins, config);
        
        % calculate the maximum map size. for the randomized dimension, the
        % number of bins for this dimension is found using the first set of
        % argins.randomBinEdges.
        mapSize(i) = length(argins.randomBinEdges{1}) - 1;
    end
end

dimInds = [dimInds, data_behav.Trial]; % append trial IDs to dimInds

% dimension indices outside of the range of binEdges were set to NaNs.
% find which frames in dimInds are NaNs and remove them from the data.
nanInd = sum(isnan(dimInds), 2) > 0; % NaN index for frames (or rows)
data_ROIs(nanInd, :) = [];  % delete nanInd rows in ROI DFF
data_behav(nanInd, :) = []; % delete nanInd rows in behavioral data
dimInds(nanInd, :) = [];    % delete nanInd rows in dimInds

clear nanInd vals* vAllTrials i
wb = waitbar(1, wb, 'Data prepared.');
%% REAL DATA: calculate the lambda(x), p(x), and MI for each ROI

wb = waitbar(0, wb, 'Beginning real data analysis.');
if argins.calcRealData == 1 % do the real data (not shuffled)
    unfilteredAverageMaps(numROIs).averageMap = []; % preallocate field
    % if a dimension is 'Evidence', then get the range of evidence values
    % from the behavioral data. this is useful for plotting. this does not
    % necessarily agree with what the user specifies for argins.binEdges.
    if any(strcmp(argins.dimensions, 'Evidence'))
        % get the min and max evidence value
        unfilteredAverageMaps(1).evidenceRange = ...
            [min(data_behav.Evidence), max(data_behav.Evidence)];
    end
    
    % preallocate p(x) using mapSize. because the analysis is written for 1
    % to 3 dimensions, the code has different logic with similar operations
    % to make it compatible with the number of dimensions.
    if numDimensions == 1 % 1D pX
        pX = zeros(mapSize, 1);
    else % 2D or 3D pX
        pX = zeros(mapSize);
    end
    
    % calculate lambda(x) for each ROI
    for ROI = 1:numROIs % loop through each ROI
        v = data_ROIs(:, ROI); % get the DFF values for an ROI
        
        % do log transformation, if applicable
        if strcmp(argins.log, 'log10')
            v = log10(v + 1);
        end
        
        % calculate lambda(x). pX is the same for every ROI because it
        % describes the number of frames that the animal is in each bin
        % across the data, which doesn *not* vary across ROIs. pX is
        % calculated along with lambda(x) from helper_getSkaggs_createMaps
        [unfilteredAverageMaps, pX] = helper_getSkaggs_createMaps(v, ...
            dimInds, ROI, numDimensions, unfilteredAverageMaps, pX, ...
            mapSize);
        
        % update waitbar
        wb = waitbar(ROI/numROIs, wb, 'Constructing 2D Maps for Real Data');
    end
    clear v allDimensions
    
    % calculate MI for each ROI
    skaggsMetric.skaggs_real = zeros(1, numROIs); % MI metric for each ROI
    skaggsMetric.rMeans_real = zeros(1, numROIs); % mean DFF per bin for each ROI
    
    for ROI = 1:numROIs % loop through ROIs
        % calculate MI metric and rMeans for each ROI
        [skaggsMetric.skaggs_real(ROI), skaggsMetric.rMeans_real(ROI)] =...
            helper_getSkaggs_calcSkaggs(unfilteredAverageMaps, ROI, ...
            argins, config, pX);
    end
    % update waitbar
    wb = waitbar(1, wb, 'Real data finished.');
end

%% SHUFFLE DATA: retrieve/write shuffle tests

% for the shuffle data, each shuffle test is saved to the path indicated in
% input6shuffleFilename. This saves time on subsequent calls of the
% function in the future. When the real data is compared to the shuffle
% data, the argins structures are compared to ensure that the
% parameters between the two are the same.

if argins.numberShuffles > 0
    rng(config.shuffleRNG); % set RNG for reproducibility
    
    % shuffleAvgMapSize contains info about number of bins, shuffles, ROIs
    shuffleAvgMapSize = [mapSize, argins.numberShuffles, numROIs];
    
    % preallocate shuffleAvgMaps that contains a lambda(x) for each shuffle
    % test for an ROI.
    shuffleAvgMaps = zeros(shuffleAvgMapSize);
    
    % preallocate fields for shuffle MI metrics
    skaggsMetric.skaggs_shuffles = NaN(argins.numberShuffles, numROIs);
    skaggsMetric.rMeans_shuffles = NaN(argins.numberShuffles, numROIs);
    
    % check if a shuffle test has already been written to folder. if so,
    % load it. if not, create a new shuffle.
    for s = 1:argins.numberShuffles % loop through shuffles
        
        shuffleFound = []; % preallocate shuffle results
        wb = waitbar(s/argins.numberShuffles, wb, ...
            'Checking for Shuffle Tests'); % update waitbar
        try
            % try to load a shuffle using the shuffleFilename
            shuffleFound = load(sprintf(argins.shuffleFilename, s));
            % assign shuffle test to shuffleFound
            shuffleFound = shuffleFound.shuffleFound;
        catch % if no shuffle is found, then continue.
        end
        
        % if a shuffle has been loaded, check the shuffle argins to
        % ensure it agrees with the current function call
        if ~isempty(shuffleFound) % a shuffle has been loaded
            % if the argins do not agree, throw an error.
            if ~isequal(shuffleFound.argins, argins)
                disp('Shuffled Input Arguments:');
                disp(shuffleFound.argins);      % print shuffle argins
                disp('Current Input Arguments:');
                disp(argins);                   % print current argins
                error(['Check saved shuffle. Shuffled input arguments', ...
                    ' do not match the current input arguments.']);
            end
        end
        
        %% if no shuffle was loaded, then write a new shuffle test.
        if isempty(shuffleFound) % no shuffle was loaded
            % update waitbar
            wb = waitbar(s/argins.numberShuffles, wb, 'Writing Shuffle to Folder');
            
            shuffleMani = []; % this holds the shuffled lambda(x)
            shuffleMani(numROIs).averageMap = []; %#ok<*AGROW>
            
            % for the shuffle test, do a circular shift of each ROI's DFF
            % within each trial in order to preserve relative firing rate.
            % each ROI is shuffled independently from other ROIs.
            data_ROIs_shuffle = data_ROIs; % copy DFF data
            
            for ROI = 1:numROIs % loop through ROIs
                v = data_ROIs(:, ROI); % take the DFF of one ROI
                % do log transformation, if applicable
                if strcmp(argins.log, 'log10')
                    v = log10(v + 1);
                end
                % shuffle the ROI DFF
                vShuffle = circshift_within_trials(v, dimInds(:,end), ...
                    config.circshiftIndLimit);
                % place the shuffled DFF back into
                data_ROIs_shuffle(:, ROI) = vShuffle;
                % create maps
                [shuffleMani, ~] = helper_getSkaggs_createMaps(vShuffle,...
                    dimInds, ROI, numDimensions, shuffleMani, pX, mapSize);
            end
            
            % save shuffled results in the structure, shuffleFound
            shuffleFound.shuffleMani = shuffleMani; % shuffled lambda(x)
            shuffleFound.argins = argins; % shuffled argins
            shuffleFound.shuffledDFF = data_ROIs_shuffle; % shuffled DFF
            shuffleFound.config = config; % shuffled config
            
            % create filename for shuffled data file
            filename = sprintf(argins.shuffleFilename, s);
            % save shuffled data file
            save(filename, 'shuffleFound','-v7.3');
            % update waitbar
            wb = waitbar(s/argins.numberShuffles, wb, 'Writing Shuffles: Saved');
        end
        
        %% calculate MI metric for shuffle data
        if argins.calcRealData == 1
            % because we need the real data's pX, we only calculate shuffle
            % skaggs metric when we calculate the real data
            for ROI = 1:numROIs % loop through ROIs
                % find MI value from shuffleFound.shuffleMani
                [skaggsMetric.skaggs_shuffles(s, ROI), ...
                    skaggsMetric.rMeans_shuffles(s, ROI), ~, rSmooth] = ...
                    helper_getSkaggs_calcSkaggs(...
                    shuffleFound.shuffleMani, ROI, shuffleFound.argins, ...
                    shuffleFound.config, pX);
                
                % keep the rSmooth in shuffleAvgMaps used for the shuffle
                % vs real lambda(x) comparisons (pixelwise).
                switch numDimensions
                    case 1 % 1D case
                        shuffleAvgMaps(:, s, ROI) = rSmooth;
                    case 2 % 2D case
                        shuffleAvgMaps(:, :, s, ROI) = rSmooth;
                    case 3 % 3D case
                        shuffleAvgMaps(:, :, :, s, ROI) = rSmooth;
                end
            end
        end
        
        % update wait bar
        wb = waitbar(s/argins.numberShuffles, wb, ...
            'Mutual Information Metric Computed.');
    end
    
    %% compare the real MI values to the shuffle distribution
    % take the mean of the MI values across shuffles
    skaggs_shuffles_mean = nanmean(skaggsMetric.skaggs_shuffles);
    % take the standard deviation across shuffles
    skaggs_shuffles_std = nanstd(skaggsMetric.skaggs_shuffles);
    
    % preallocate list to hold ROIs with significantly higher MI than
    % shuffles
    skaggsSignificantROIs = [];
    for ROI = 1:numROIs % loop through ROIs
        % if the real MI value for an ROI is greater than shuffle mean + N
        % standard deviations where N is argins.significanceThreshold,
        % then that ROI passes shuffle test
        if (skaggsMetric.skaggs_real(ROI)) > (skaggs_shuffles_mean(ROI) ...
                + argins.significanceThreshold*skaggs_shuffles_std(ROI))
            % save ROI as significant
            skaggsSignificantROIs = [skaggsSignificantROIs, ROI];
        end
    end
    % place significant MI ROIs in output field
    skaggsMetric.sigROIs = skaggsSignificantROIs;
end

clear rSmooth shuffleAvgMapSize skaggs_shuffles*


%% Pixelwise: ROI receptive field computation from lambda(x)
% in this analysis, the lambda(x) for the real data is compared to the
% lambda(x) for the shuffle data. with N shuffles, each bin for x has a
% distribution of N shuffle test values where the mean and standard
% deviation are found.

% do this analysis if the real and shuffle analyses were both done
if argins.calcRealData == 1 && argins.numberShuffles > 0
    
    % preallocate fields
    pixelwise(numROIs).avgRealMap = []; % real lambda(x)
    pixelwise(numROIs).avgShuffleMap = []; % mean lambda(x) across shuffles
    pixelwise(numROIs).stdShuffleMap = []; % std lambda(x) across shuffles
    pixelwise(numROIs).thresholdMap = []; % mean shuffle lambda(x) + (config.shroudThreshold * std lambda(x))
    pixelwise(numROIs).realPassThreshold = []; % significant bins where real lambda(x) is over thresholdMap
    pixelwise(numROIs).realSignifMap = []; % avgRealMap only including significant bins
    
    for ROI = 1:numROIs % loop through ROIs
        % do the bin-wise comparison of the real lambda(x) with shuffled
        % lambda(x)
        pixelwise = helper_getSkaggs_pixelwise(unfilteredAverageMaps, ...
            argins, config, ROI, pixelwise, numDimensions, shuffleAvgMaps);
        % update waitbar
        wb = waitbar(ROI/numROIs, wb, 'Creating Pixelwise');
    end
    clear ROI s n shuffleAvgMaps
end

%% Compare each shuffle test to the rest of the shuffle distribution
% in order to find the chance number of MI significant ROIs in shuffles.
% for each shuffle, find the number of MI significant ROIs by treating
% that shuffle as the real data. for example, for a distribution of 50
% shuffles, take shuffle test 1 and find the number of MI significant ROIs
% when compared to the distribution of shuffles 2-50. then take shuffle
% test 2 and find the number of MI significant ROIs when compared to the
% distribution of shuffles 1 and 3-50. each shuffle test will then have a
% number of MI significant ROIs found by chance, and this distribution from
% the shuffle data is then compared to the number of MI significant ROIs in
% the real data.

if argins.calcRealData == 1 % only applicable if real data is calculated
    if argins.numberShuffles > 0 % only applicable if shuffles were done
        % preallocate shuffleComparison. holds number of MI significant
        % ROIs in the shuffle data
        shuffleComparisons = NaN(argins.numberShuffles, 1);
        for s = 1:argins.numberShuffles % loop through shuffles
            nontarget = 1:argins.numberShuffles; % shuffle distribution
            nontarget(s) = []; % take shuffle 's' out of distribution
            
            % recalculate the distribution's mean and std
            nontarget_avg = mean(skaggsMetric.skaggs_shuffles(nontarget, :));
            nontarget_std = std(skaggsMetric.skaggs_shuffles(nontarget, :));
            
            % find the significance threshold
            nontarget_thresh = nontarget_avg + ...
                argins.significanceThreshold * nontarget_std;
            
            % save the number of significant ROIs found by chance for this
            % shuffle 's'
            shuffleComparisons(s) = sum(...
                skaggsMetric.skaggs_shuffles(s, :) > nontarget_thresh);
        end
        % save shuffleComparisons in output
        skaggsMetric.shuffleComparisons_sigROIs = shuffleComparisons;
        
        % compare number of MI significant ROIs from real data to shuffles
        shuffleComparisons_mean = mean(shuffleComparisons); % shuffle mean
        shuffleComparisons_std = std(shuffleComparisons); % shuffle std
        
        % if significantly more MI ROIs in real data than shuffle, then the
        % real MI ROIs are above the shuffle distribution
        if length(skaggsSignificantROIs) > shuffleComparisons_mean + ...
                argins.significanceThreshold * shuffleComparisons_std
            skaggsMetric.sigROIsAboveShuffleDistribution = 1;
        else % if not, the real data is not above the shuffle distribution
            skaggsMetric.sigROIsAboveShuffleDistribution = 0;
        end
    end
end
clear shuffleComparisons* nontarget* m s

%% Create sequence plots
% sort the ROIs' lambda(x) based on their peak mean DFF to see if a
% sequence is made where the ROIs' activities tile the binned space. a
% cross-validation procedure is also used whereby for a given ROI, its
% lambda(x) is recalculated for all trials split 50/50 into a training and
% testing set. the sequence is sorted for all ROIs on the training set, and
% this same sorting index is applied to the testing set to see if the
% sequence is conserved. the sequence plot analysis is only valid for
% 1D analyses. in 2D or 3D analyses, getSkaggs() is called individually for
% each dimension.

if argins.calcRealData == 1 % only applicable if real data is analyzed
    
    % if no shuffles were used, then no MI significance test was done.
    % instead of using MI significant ROIs, use all ROIs
    if argins.numberShuffles == 0
        skaggsSignificantROIs = 1:numROIs; % use all ROIs
    end
    
    % if numDimensions == 1 and the dimension is not binary, find the
    % sequence plot. a sequence for a binary variable that only contains
    % two possible values doesn't seem interpretable
    if numDimensions == 1 && ... % 1D case and no binary variables
            any(strcmp(argins.dimensions{1}, {'Choice', ...
            'ChoiceCorrect', 'PriorChoice', 'PriorCorrect', ...
            'NearestCueType', 'PrevEvidence'})) == 0
        
        % calculate the sequence plot and save outputs
        [~, sequencePlot(1).CV.matrixTrain, ...      % training plot
            sequencePlot(1).CV.matrixTest, ...       % testing plot
            sequencePlot(1).noCV.unsortedMatrix, ... % unsorted no CV plot
            sequencePlot(1).noCV.sortedMatrix, ...   % sorted CV plot
            sortInd_odd, ...                      % training set sort index
            ~, sortInd_noCV, ...                  % no CV set sort index
            sequencePlot(1).shuffleTest.realCorr, ...% CV kendalls corr
            sequencePlot(1).ROIsPassedCV] ...        % ROIs that passed CV
            = helper_getSkaggs_CVseqplot(...
            dimInds, skaggsSignificantROIs, instTrialsAverageMap, ...
            data_ROIs, pX, argins, config);
        
        % save dimension names
        sequencePlot(1).CV.dimensionTrain = argins.dimensions{1};
        sequencePlot(1).CV.dimensionTest = argins.dimensions{1};
        sequencePlot(1).noCV.dimension = argins.dimensions{1};
        
        % save ROI labels of sorted plots
        sequencePlot(1).CV.sortedROIlabelsTrain = ... % CV train set
            sequencePlot(1).ROIsPassedCV(sortInd_odd);
        
        sequencePlot(1).CV.sortedROIlabelsTest = ... % CV test set
            sequencePlot(1).ROIsPassedCV(sortInd_odd);
        
        sequencePlot(1).noCV.sortedROIlabels = ... % no CV set
            sequencePlot(1).ROIsPassedCV(sortInd_noCV);
        
        % create sequence plots for each shuffle
        if argins.numberShuffles > 0
            % preallocate variable to hold correlation values
            kendallTauShuffle = NaN(argins.numberShuffles, 1);
            for i = 1:argins.numberShuffles % loop through shuffles
                % load shuffle data
                load(sprintf(argins.shuffleFilename, i), 'shuffleFound');
                
                % get shuffled ROI DFF
                data_ROIs_shuffle = shuffleFound.shuffledDFF;
                
                % only compute the kendalls correlation for shuffle data
                [~, ~, ~, ~, ~, ~, ~, ~, kendallTauShuffle(i)] ...
                    = helper_getSkaggs_CVseqplot(dimInds, ...
                    skaggsSignificantROIs, instTrialsAverageMap, ...
                    data_ROIs_shuffle, pX, argins, config);
                
                % update waitbar
                wb = waitbar(i/argins.numberShuffles, wb, ...
                    'Comparing Sequence Plot to Shuffles');
            end
            % save shuffled kendalls correlation
            sequencePlot(1).shuffleTest.shuffledCorr = kendallTauShuffle;
        end
    end
    clear shuffleFound sortInd_* data_ROIs_shuffle kendallTauShuffle
    clear instTrialsAverageMap
    
    % Create sequence plots for 2D or 3D analyses. call getSkaggs() again
    % for each dimension using the MI significant ROIs
    if numDimensions > 1 % 2D or 3D case
        for i = 1:numDimensions % loop through each dimension
            % call getSkaggs for each dimension using many of the same
            % input arguments as the original function call. the following
            % arguments change:
            %
            % input5numberShuffles - change to 0. no shuffle tests done
            % input6shuffleFilename - changed so it does not conflict with
            % current function call's shuffle filename. because no shuffles
            % are made, the name isn't used for anything but still has to
            % end in '%d.mat'
            % input15whichROIs - changed to MI significant ROIs only
            % instead of the ROIs originally specified.
            recurOut = getSkaggs(...
                argins.DFFthresholding, argins.log, ...
                argins.gaussianSmoothingSigma, argins.trialsAverageMap, ...
                0, 'unused_%d.mat', 1, argins.imagingFilename, ...
                argins.significanceThreshold, argins.dimensions(i), ...
                argins.binEdges(i), argins.randomDimMethod, ...
                argins.randomBinEdges, argins.taskType, ...
                skaggsSignificantROIs, argins.gaussianSmoothingToggle);
            
            % the sequence plot analysis does not apply for binary
            % variables so skip those binary dimensions
            try
                % each dimension has its own sequence plot. copy the
                % relevant information into the output structure.
                
                % dimension name
                sequencePlot(i).dimension = argins.dimensions{i};
                % CV sequence plot
                sequencePlot(i).CV = recurOut.sequencePlot.CV;
                % no CV sequence plot
                sequencePlot(i).noCV = recurOut.sequencePlot.noCV;
                
                % copy the ROI sorting labels. because these ROI labels
                % were found in recurOut, they are redescribed in
                % terms of the original MI significant ROIs
                sequencePlot(i).CV.sortedROIlabelsTrain = ... % CV train
                    skaggsSignificantROIs(...
                    sequencePlot(i).CV.sortedROIlabelsTrain);
                sequencePlot(i).CV.sortedROIlabelsTest = ... % CV test
                    skaggsSignificantROIs(...
                    sequencePlot(i).CV.sortedROIlabelsTest);
                sequencePlot(i).noCV.sortedROIlabels = ... % no CV
                    skaggsSignificantROIs(...
                    sequencePlot(i).noCV.sortedROIlabels);
                sequencePlot(i).ROIsPassedCV = ... % ROIs passed CV
                    skaggsSignificantROIs(...
                    recurOut.sequencePlot.ROIsPassedCV);
            catch
            end
        end
    end
end

out.unfilteredAverageMaps = unfilteredAverageMaps; % ROIs' lambda(x)
out.pixelwise = pixelwise; % real vs shuffled lambda(x)
out.pX = pX; % p(x), the number of frames the animal spent in each bin
out.skaggsMetric = skaggsMetric; % MI analysis
out.sequencePlot = sequencePlot; % sequence plot analysis
out.config = config; % configurations structure
out.argins = argins; % input arguments structure
out.DFF = data_ROIs; % ROI DFF data
out.behav = data_behav; % behavioral data
out.dimensionIndices = dimInds(:, 1:numDimensions); % bin index per frame
out.ROIs = 1:size(data_ROIs,2); % number of ROIs

delete(wb);
end

%% ******* Helper functions ******* %%
%% CHECK INPUTS
function [] = helper_getSkaggs_checkInputs(argins)
%
% Sample call: helper_getSkaggs_checkInputs(argins);
%
% Check several input arguments into getSkaggs() that may not be specified
% correctly. If an argument is incorrect, then an error is thrown. This
% function only checks common mistakes, not all possible errors from bad
% arguments input into getSkaggs().
%
% ***** INPUTS *****
% argins - structure of input arguments into getSkaggs()


if length(argins.DFFthresholding) ~= 2
    % input1 needs an array of length two for DFF thresholding.
    % first element controls the window size.
    % second element controls thresholding for data (x * robustSTD)
    error(['Input 1 (DFF Thresholding) must be specified as an array', ...
        ' of length 2.']);
end
if sum(strcmp(argins.log, {'log10', 'noLog'})) ~= 1
    % input2 needs either 'log10' or 'noLog'.
    % if user specied neither or both, throw an error.
    error('Input 2 (Log) must be specified as "noLog" or "log10".');
end
if argins.gaussianSmoothingSigma < 0
    % input3 sigma for smoothing with a Gaussian filter
    % must be greater than 0.
    error('Input 3 (Gaussian Smoothing) must be greater than 0.');
end
if (floor(argins.numberShuffles) ~= argins.numberShuffles) || ...
        argins.numberShuffles < 0
    % input5 must be an integer and greater than zero
    error('Input 5 (Number of Shuffles) must be a positive integer.');
end
if ~contains(argins.shuffleFilename, '%d.mat')
    % input6 must end in '%d.mat' to save shuffles on each iteration of
    % shuffle for-loop
    error('Input 6 (Shuffle Filename) must end with "%d.mat".');
end
if sum(argins.calcRealData ~= [0 1]) ~= 1
    % calculate real data MI values must be 0 (don't calculate it) or 1
    % (calculate it).
    error('Input 7 (Calculate Real Data) must be 0 or 1.');
end
if argins.significanceThreshold < 0
    % this is the threshold for the MI significance test, so it must be
    % greater than zero.
    error('Input 9 (Significance Threshold) must be greater than 0.');
end
if ~iscell(argins.dimensions) || ~iscell(argins.binEdges)
    % for consistency among 1D, 2D, and 3D MI analyses, input10 and input11
    % must be specified in cell arrays.
    error(['Input 10 (Dimensions) and Input 11 (Bin Edges) ', ...
        'must be specified in cell arrays.']);
end
if length(argins.dimensions) ~= length(argins.binEdges)
    % for each dimension, a set of bin edges is needed.
    % if the lengths of the cell arrays aren't equal, throw an error.
    error(['Input 10 (Dimensions) and Input 11 (Bin Edges) ', ...
        'must be the same length.']);
end
if (length(argins.randomDimMethod) > 3 || ...
        ~iscell(argins.randomDimMethod)) || ...
        (length(argins.randomBinEdges) > 3 || ...
        ~iscell(argins.randomBinEdges))
    % because the random dimension method samples dimension 1 from the
    % joint distribution with dimension 2, dimension 1 and 2 must be
    % specified, each with a set of bin edges.
    % input12 and input13 need to be cell arrays of length 2
    
    % EDIT: 2021/1/14 - changed to do 3D random dimension too, with some
    % functionality, like the plotting, removed
    error(['Input 12 (Random Dim Method) and Input 13 (Random Bin ', ...
        'Edges) must be a cell array less than length 3.']);
end
if sum(strcmp(argins.taskType, {'towers', 'alternationJeff'})) ~= 1
    % the datafiles are either the accumulating towers task or the
    % alternation task
    error(['Input 14 (Task Type) must be specified as either ', ...
        'towers or alternation.']);
end
if sum(strcmp(argins.gaussianSmoothingToggle, ...
        {'both', 'skaggsOnly', 'pixelwiseOnly'})) ~= 1
    % gaussian smoothing may occur for the MI analysis and/or the pixelwise
    % analysis.
    error(['Input16 (Gaussian Smoothing Toggle must be either ', ...
        '"both", "skaggsOnly", or "pixelwiseOnly".']);
end

% check dimensions in the analysis. if all specified dimensions do *not*
% vary in a trial (Choice, PriorChoice, ChoiceCorrect, PriorCorrect,
% TrialProb), then the shuffling method (circshift_within_trial) is *not*
% appropriate. These variables would *not* be shuffled after the shuffling
% method. If there is at least one variable that varies in each trial, then
% the shuffling test is appropriate because it shuffles the joint
% distribution of these variables.

c = 0; % save the number of dimensions that do *not* vary in a trial
for i = 1:length(argins.dimensions) % check all dimensions
    if any(strcmp(argins.dimensions{i}, {'Choice', 'PriorCorrect', ...
            'PriorChoice', 'ChoiceCorrect', 'TrialProb'}))
        c = c+1; % increase c if the dimension doesn't vary in a trial
    end
end
% check if there's at least one dimension that *does* vary within a trial
% if not, then throw an error.
if c == length(argins.dimensions)
    error(['The input dimension(s) is not compatible with the ', ...
        'shuffling method (circular shift within trials).']);
end


end

%% DISCRETIZE INDICES
function out = helper_getSkaggs_discretize(vals, argins, config, i, randomIndex)
%
% Sample Call: vals2 = helper_getSkaggs_discretize(vals, argins, config, i, 0);
%
% discretize vals based on the bin edges in argins.
%
% ***** INPUTS *****
%
% vals - 1-by-N vector of values to discretize
% argins - input arguments structure that contains 'binEdges' field
% config - configuration structure that contains default number of bins
% i - index specifying the i-th dimension
% randomIndex - double specifying whether the dimension is non-randomized
%   (randomIndex == 0) or randomized (randomIndex == 1 or 2).
%
% ***** OUTPUTS *****
%
% out - 1-by-N vector of discretized input 'vals'
%

out = []; % preallocate output
if randomIndex == 0 % if randomIndex == 0, then it's not randomized
    if ~isempty(argins.binEdges{i}) % bin edges are specified
        out = discretize(vals, argins.binEdges{i}); % discretize using bin edges
    else % bin edges are *not* specified, discretize using default number of bins
        out = discretize(vals, config.discretizeDefaultNumberBins);
    end
elseif any(randomIndex == [1, 2, 3])
    % if randomIndex == 1 or 2, then it's a randomized dimension 1 (randomIndex == 1)
    % sampled with respect to dimension 2 (randomIndex == 2)
    if ~isempty(argins.randomBinEdges{randomIndex}) % bin edges are specified
        out = discretize(vals, argins.randomBinEdges{randomIndex});
    else % bin edges are *not* specified, discretize using default number of bins
        out = discretize(vals, config.discretizeDefaultNumberBins);
    end
end
end

%% RANDOM DIMENSION
function binnedR = helper_getSkaggs_randomDimension(binnedDim1, binnedDim2)
%
% Sample Call: binnedR = helper_getSkaggs_randomDimension(binnedDim1, binnedDim2);
%
% Sample uniformly from the joint distribution of dimension 1 w.r.t. dimension 2
%
% ***** INPUTS *****
%
% binnedDim1 - 1-by-N vector of binned values for dimension 1
% binnedDim2 - 1-by-N vector of binned values for dimension 2
%
% ***** OUTPUTS *****
%
% binnedR - 1-by-N vector of dimension 1 values sampled w.r.t. dimension 2
%

binnedR = NaN(size(binnedDim2)); % preallocate output
% for each bin of dim 2, find the corresponding values of dim1 and sample
% from those values.
for i = min(binnedDim2):max(binnedDim2) % loop through dim 2 bins
    jointInd = binnedDim2 == i; % index where dim2 is equal to loop index
    % sample from the binnedDim1 values where binnedDim2 is equal to i.
    jointDist = datasample(binnedDim1(jointInd), ...
        length(binnedDim1(jointInd)));
    % place sampled observations into binnedR where jointInd is true.
    binnedR(jointInd) = jointDist;
end
end

%% RANDOM DIMENSION (for 3 Dimensions)
function binnedR = helper_getSkaggs_randomDimension3D(binnedDim1, binnedDim2, binnedDim3)
%
% Sample Call: binnedR = helper_getSkaggs_randomDimension(binnedDim1, binnedDim2, binnedDim3);
%
% Sample uniformly from the joint distribution of dimension 1 w.r.t.
% dimensions 2 and 3
%
% ***** INPUTS *****
%
% binnedDim1 - 1-by-N vector of binned values for dimension 1
% binnedDim2 - 1-by-N vector of binned values for dimension 2
% binnedDim3 - 1-by-N vector of binned values for dimension 3
%
% ***** OUTPUTS *****
%
% binnedR - 1-by-N vector of dimension 1 values sampled w.r.t. dimensions 2
% and 3
%

binnedR = NaN(size(binnedDim2)); % preallocate output
% for each bin of dim 2, find the corresponding values of dim1 and sample
% from those values.
for i = min(binnedDim2):max(binnedDim2) % loop through dim 2 bins
    for j = min(binnedDim3):max(binnedDim3)
        
        jointInd = binnedDim2 == i & binnedDim3 ==j; % index where dim2 and dim3 are equal to loop indices
        % sample from the binnedDim1 values where binnedDim2 is equal to i and binnedDim3 is equal to j.
        jointDist = datasample(binnedDim1(jointInd), ...
            length(binnedDim1(jointInd)));
        % place sampled observations into binnedR where jointInd is true.
        binnedR(jointInd) = jointDist;
    end
end
end

%% PLOT RANDOM VARIABLES
function [] = plot_randomVariables_pX(binnedDim1, binnedDim2, ...
    binnedDim1_R, argins, config)
% creates heatmap of pX of randomized dimensions compared to real
% dimensions
%
% ***** INPUTS *****
% binnedDim1 - vector of indices for real dimension 1
% binnedDim2 - vector of indices for real dimension 2
% binnedDim1_R - vector of indices for randomized dimension 1
% binnedDim2_R - vector of indices for randomized dimension 2
% argins - input arguments structure containing names for these dimensions
%
% ***** OUTPUTS ******
% 1-by-2 subplot figure. the first plot shows the pX of the non-randomized
% dimensions. the second plot shows the pX of the randomized dimensions.
%

nBins1 = length(argins.randomBinEdges{1}) - 1; % number of bins for dim1
nBins2 = length(argins.randomBinEdges{2}) - 1; % number of bins for dim2

figure; % create new figure

%% first plot: pX of non-randomized dimensions
subplot(1,2,1);

% use helper_getSkaggs_createMaps to get pX for non-randomized dimension 1
% w.r.t. dimension 2.
v = zeros(size(binnedDim1)); % required as input to helper_getSkaggs_createMaps but unused in plot_randomVariables_pX
mDims = [binnedDim1, binnedDim2]; % dimension indices
inst(1).averageMap = []; % required as input to helper_getSkaggs_createMaps but unused in plot_randomVariables_pX
pX_nonRandomized = zeros(nBins1, nBins2); % sets size of pX

% calculate pX
[~, pX_nonRandomized] = helper_getSkaggs_createMaps(v, mDims, 1, 2, ...
    inst, pX_nonRandomized, [nBins1, nBins2]);

% plot pX
imagesc(pX_nonRandomized');

% label title
title('Non-randomized p(x)');

%% second plot: pX of randomized dimensions

subplot(1,2,2);
% use helper_getSkaggs_createMaps to get pX for *randomized* dimension 1
% w.r.t. dimension 2.
mDims = [binnedDim1_R, binnedDim2]; % dimension indices for randomized dim
pX_randomized = zeros(nBins1, nBins2); % sets size of pX

% calculate pX
[~, pX_randomized] = helper_getSkaggs_createMaps(v, mDims, 1, 2, ...
    inst, pX_randomized, [nBins1, nBins2]);

% plot pX
imagesc(pX_randomized');

% label title
title('Randomized p(x)');

%% add axis labels and colorbar
for i = 1:2 % loop through each subplot
    subplot(1,2,i);
    % label axes
    xlabel(argins.randomDimMethod(1));
    ylabel(argins.randomDimMethod(2));
    % for visualization, set colorbar range from 0 to config.random_pX_caxisLim
    caxis([0, config.random_pX_caxisLim]);
    colorbar;
    % change tick mark orientation
    xtickangle(90);
end

end

%% CREATE MAPS
function [sMaps, pXinst] = helper_getSkaggs_createMaps(v, dimInds, ROI, ...
    numDimensions, sMaps, pX, mapSize)
% fills sMaps structure with lambda(x), the mean DFF across bins, for ROIs
%
% ***** INPUTS *****
% v - vector of DFF values for an ROI
% dimInds - matrix of dimension indices
% ROI - index for which ROI to do
% numDimensions - double specifying the number of dimensions
% sMaps - structure to hold the mean DFF across bins for ROIs
% pX - matrix showing the number of frames the animal spent in each bin
% mapSize - the number of bins for each dimension in pX and sMaps fields
%
% ***** OUTPUTS ******
% sMaps - 1-by-numROIs structure that holds lambda(x) for each ROI.
% pXinst - pX, the number of frames that the animal spends in each bin.

pXinst = zeros(size(pX)); % preallocate pXinst

%% 1D lambda(x) and p(x)
if numDimensions == 1
    sMaps(ROI).averageMap = NaN([mapSize, 1]); % preallocate lambda(x)
    for i = min(dimInds(:,1)):max(dimInds(:,1)) % first dimension loop
        vInd = dimInds(:, 1) == i; % find where dimInds equals 1D bin index
        sMaps(ROI).averageMap(i) = mean(v(vInd)); % get mean DFF for bin
        pXinst(i) = sum(vInd); % pXinst is the number of frames equal to the bin index
    end
end

%% 2D lambda(x) and p(x)
if numDimensions == 2
    sMaps(ROI).averageMap = NaN(mapSize); % preallocate lambda(x)
    for i = min(dimInds(:,1)):max(dimInds(:,1)) % first dimension loop
        for j = min(dimInds(:,2)):max(dimInds(:,2)) % second dimension loop
            vInd = dimInds(:, 1) == i & dimInds(:,2) == j; % find where dimInds equals 2D bin index
            sMaps(ROI).averageMap(i, j) = mean(v(vInd)); % get mean DFF for bin
            pXinst(i, j) = sum(vInd); % pXinst is the number of frames equal to the bin index
        end
    end
end

%% 3D lambda(x) and p(x)
if numDimensions == 3
    sMaps(ROI).averageMap = NaN(mapSize); % preallocate lambda(x)
    for i = min(dimInds(:,1)):max(dimInds(:,1)) % first dimension loop
        for j = min(dimInds(:,2)):max(dimInds(:,2)) % second dimension loop
            for k = min(dimInds(:,3)):max(dimInds(:,3)) % third dimension loop
                vInd = (dimInds(:, 1) == i & dimInds(:,2) == j) & dimInds(:,3) == k; % find where dimInds equals 3D bin index
                sMaps(ROI).averageMap(i, j, k) = mean(v(vInd)); % get mean DFF for bin
                pXinst(i, j, k) = sum(vInd); % pXinst is the number of frames equal to the bin index
            end
        end
    end
end

end

%% SMOOTHING LAMBDA(X)
function [rSmooth] = helper_getSkaggs_smoothing(r, argins, config)
%
% Smooths r by convolution with a Gaussian filter.
%
% if the sigma for Gaussian smoothing (argins.gaussianSmoothingSigma) is
% greater than 0, then smooth lambda(x) before finding
% the Mutual Information metric. the smoothing method is a convolution of
% lambda(x) with a Gaussian filter (sigma is specified by
% argins.gaussianSmoothingSigma, while the size of the smoothing window is
% set by config.gaussianSmoothingParameters). the smoothing method is only
% applicable for 1D and 2D analyses. NaNs in the lambda(x) are preserved in
% the smoothed version by using nanconv().
%
% ***** INPUTS *****
% r - lambda(x) to be smoothed
% argins - structure containing input arguments from getSkaggs()
% config - structure containing configurations
%
% ***** OUTPUTS *****
% rSmooth - smoothed r (or lambda(x))

% check if sigma is greater than zero
if argins.gaussianSmoothingSigma > 0
    switch length(argins.dimensions) % consider 1D, 2D, and 3D cases
        case 1 % 1D case
            % create 1D Gaussian filter
            imageFilter=fspecial('gaussian', ...
                [1, config.gaussianSmoothingParameters], ...
                argins.gaussianSmoothingSigma);
            
            % convolution of lambda(x) with filter.
            rSmooth = nanconv(r, imageFilter, 'edge', 'nanout', '1d');
            
        case 2 % 2D case
            % create 2D Gaussian filter
            imageFilter=fspecial('gaussian', ...
                [config.gaussianSmoothingParameters, ...
                config.gaussianSmoothingParameters], ...
                argins.gaussianSmoothingSigma);
            
            % convolution of lambda(x) with filter.
            rSmooth = nanconv(r, imageFilter, 'edge', 'nanout', '2d');
        case 3 % 3D case
            if size(r,3)~=2
                % throw error. smoothing is not supported in 3D analyses.
                error('Error: 3D Gaussian Smoothing is not allowed.');
            else
                disp('BROKEN: Doing SPECIAL CASE where third dimension is 2 bins');
                
            end
            
            imageFilter=fspecial('gaussian', ...
                [config.gaussianSmoothingParameters, ...
                config.gaussianSmoothingParameters], ...
                argins.gaussianSmoothingSigma);
            
            % convolution of lambda(x) with filter.
            rSmooth(:,:,1) = nanconv(r(:,:,1), imageFilter, 'edge', 'nonanout', '2d');
            rSmooth(:,:,2) = nanconv(r(:,:,2), imageFilter, 'edge', 'nonanout', '2d');
    end
end
end

%% CALC SKAGGS
function [MImetric, rMeans, r, rSmooth] = helper_getSkaggs_calcSkaggs(...
    sMaps, ROI, argins, config, pX)
%
% calculate the Mutual Information metric (MI) for an ROI based on its
% lambda(x) and the animal's p(x)
%
% ***** INPUTS *****
% sMaps - 1-by-numROIs structure that holds lambda(x) for each ROI.
% ROI - index for which ROI to do
% argins - structure holding input arguments for getSkaggs
% config - structure holding standard configurations
% pX - matrix showing the number of frames the animal spent in each bin
%
% ***** OUTPUTS *****
%
% MImetric - mutual information value
% rMeans - mean DFF per bin
% r - smoothed (if applicable) lambda(x) without NaNs or negative values
% rSmooth - smoothed (if applicable) lambda(x) without NaNs

r = sMaps(ROI).averageMap; % get lambda(x) for an ROI

% check that smoothing toggle is either 'both' or 'skaggsOnly'
if argins.gaussianSmoothingSigma > 0 ...
        && any(strcmp(argins.gaussianSmoothingToggle, {'both', ...
        'skaggsOnly'}))
    r = helper_getSkaggs_smoothing(r, argins, config); % smooth lambda(x)
end
r(isnan(r)) = 0; % set NaNs to 0 for computeSkaggsInformation()
rSmooth = r; % save lambda(x) before thresholding

% threshold lambda(x) because computeSkaggsInformation() cannot take values
% less than zero. set negative values to zero.
r(r < 0) = 0;

% get MI metric and rMeans
[MImetric, rMeans] = computeSkaggsInformation(pX, r);
end

%% PIXELWISE
function [pixelwise] = helper_getSkaggs_pixelwise(sMaps, argins, ...
    config, ROI, pixelwise, numDimensions, shuffleAvgMaps)
% compare the lambda(x) from the real and shuffled data bin by bin.
%
% ***** INPUTS *****
% sMaps - structure containing lambda(x) for each ROI
% argins - structure containing input arguments from getSkaggs()
% config - structure containing configurations for getSkaggs()
% ROI - double, ROI label index.
% pixelwise - structure that will populated using ROI as an index
% numDimensions - double, number of dimensions
% shuffleAverageMaps - (numDimesnions+1)-D array containing lambda(x) for
% the shuffled data
%
% ***** OUTPUTS *****
% pixelwise - structure with bin-wise comparisons of lambda(x) of the real
% data with the shuffled data

r = sMaps(ROI).averageMap; % get lambda(x) for an ROI
if argins.gaussianSmoothingSigma > 0 && any(strcmp(argins.gaussianSmoothingToggle, {'both', 'pixelwiseOnly'}))
    r = helper_getSkaggs_smoothing(r, argins, config); % smooth lambda(x)
end
r(isnan(r)) = 0; % set NaNs to zero

% the field pixelwise.avgRealMap holds smoothed (if applicable)
% lambda(x) without NaNs for graphing
pixelwise(ROI).avgRealMap = r;

% for each bin in lambda(x), find the mean and standard deviation of the
% shuffle tests for that bin.
% n is the shuffle tests across bins for a given ROI.
% indexing for shuffleAvgMaps changes for different dimensions. the
% second-to-last dimension of the matrix is the shuffles, while the last
% dimension of the matrix is the ROI number
switch numDimensions % indexing changes for different dimensions
    case 1 % 1D case
        n = shuffleAvgMaps(:, :, ROI);
    case 2 % 2D case
        n = shuffleAvgMaps(:, :, :, ROI);
    case 3 % 3D case
        n = shuffleAvgMaps(:, :, :, :, ROI);
end

% using n, find the mean and standard deviation of shuffled lambda(x)
% across shuffles. compare these to the real data.

% average of shuffle lambda(x)
pixelwise(ROI).avgShuffleMap = mean(n, numDimensions + 1);
% standard deviation of shuffle lambda(x)
pixelwise(ROI).stdShuffleMap = std(n, 0, numDimensions + 1);
% significance threshold by bin. calculated using config.shroudThreshold
pixelwise(ROI).thresholdMap = pixelwise(ROI).avgShuffleMap ...
    + config.shroudThreshold*pixelwise(ROI).stdShuffleMap;
% lambda(x) bins of real data that pass significance threshold
pixelwise(ROI).realPassThreshold = pixelwise(ROI).avgRealMap ...
    > pixelwise(ROI).thresholdMap;
% real lambda(x) only containing significant bins. copy
% pixelwise.avgRealMap and set pixelwise.realPassThreshold == 0 to zero
pixelwise(ROI).realSignifMap = pixelwise(ROI).avgRealMap;
pixelwise(ROI).realSignifMap(~pixelwise(ROI).realPassThreshold) = 0;
end
%% CV SEQUENCE PLOT
function [tRank1, mTrain, mTest, mNoCVunsort, mNoCV_sort, sortInd_odd, ...
    sortInd_even, sortInd_noCV, kendallTau, ROIsPassedCV] = ...
    helper_getSkaggs_CVseqplot(dimInds, ROIs, trials, data_ROIs, pX, ...
    argins, config)
%
% Creates cross-validated and non-cross-validated sequence plots by sorting
% ROIs by the bins of their peak mean DFFs. in the CV set for each ROI,
% the trials are split into training and testing sets based on the ROI's
% maximal DFF value for a given trial, which is then converted to a rank
% after trials are sorted by their maximal DFF value. the training set
% consists of odd-ranked trials (1st, 3rd, 5th, etc. trials) whereas the
% testing set consists of even-ranked trials (2nd, 4th, 6th, etc. trials).
% only applicable for 1D analyses.
%
% ***** INPUTS *****
% dimInds - matrix of dimension indices
% ROIs - list of ROI IDs to include in sequence plot
% trials - list of trials IDs to include in sequence plot
% data_ROIs - ROI DFF data by frame
% pX - number of frames in each bin
% argins - structure containing input arguments
% config - structure containing standard configurations
%
% ***** OUTPUTS *****
% tRank1 - trial ranks based on an ROI's maximal activity in that trial
% mTrain - sorted training sequence plot
% mTest - testing sequence plot using same sorting as mTrain
% mNoCVunsort - unsorted non-CV sequence plot
% mNoCV_sort - sorted non-CV sequence plot
% sortInd_odd - sorting index used for training set to produce mTrain
% sortInd_even - sorting index if testing set was sorted independently of
%                   the training set
% sortInd_noCV - sorting index used for noCV set.
%                 mNoCVunsort(sortInd_noCV, :) gives mNoCV_sort
% kendallTau - Kendall's correlation coefficient between sortInd_odd and
%               sortInd_even
% ROIsPassedCV - list of ROIs that had significant correlation values (p <
%                   0.05) between training and testing sequence plots

ind = dimInds(:, 1);        % bin indices
vTrials = dimInds(:, end);  % trial indices

% preallocate sequence plots and tRank1
%% change made here on 6/22
% mOdd = NaN(length(ROIs), length(min(ind):max(ind)));    % training plot
mOdd = NaN(length(ROIs), length(pX));    % training plot
%%
mEven = NaN(size(mOdd));                                % testing plot
mNoCV = NaN(size(mOdd));                                % no CV plot
tRank1 = zeros(length(ROIs), length(trials));           % trial ranks

for j = 1:length(ROIs) % loop through ROIs
    % for initial sorting, only look at odd ranks (1st, 3rd, 5th, etc)
    % as the training set. trialDFF is the criterion value for the trials
    % used to determine a trial's rank
    trialDFF = zeros(length(trials), 1); % preallocate trial DFF
    
    for t = 1:length(trialDFF) % loop through trials
        % get ROI DFF values for that trial
        vals = data_ROIs(ismember(vTrials, trials(t)), ROIs(j));
        % save the maximum DFF value for that trial
        trialDFF(t) = nanmax(vals);
    end
    
    [~, tInd] = sort(trialDFF, 'descend'); % sort the trials by max DFF
    tRank1(j, tInd) = 1:length(trialDFF); % convert sort indices into ranks
    
    % ROI DFF activity. here, 'odd' means the odd-ranked trials that are
    % ordered in tInd. 'even' means the even-ranked trials in tInd.
    % similar computations will occur for training, testing, and no CV sets
    vOdd  = data_ROIs(ismember(vTrials, trials(tInd(1:2:end))), ROIs(j)); % training set DFF
    vEven = data_ROIs(ismember(vTrials, trials(tInd(2:2:end))), ROIs(j)); % testing set DFF
    vNoCV = data_ROIs(ismember(vTrials, trials), ROIs(j)); % no CV set DFF
    
    % bin indices associated with ROI DFF frames
    vInd_odd  = ind(ismember(vTrials, trials(tInd(1:2:end)))); % training set bin indices
    vInd_even = ind(ismember(vTrials, trials(tInd(2:2:end)))); % testing set bin indices
    vInd_noCV = ind(ismember(vTrials, trials)); %  no CV set bin indices
    
    % for each bin value, find the ROI's mean DFF for that bin in the
    % training, testing, and no CV sets
    for k = min(ind):max(ind) % loop through bin values
        mOdd(j, k) = mean(vOdd(vInd_odd == k)); % training set
        mEven(j,k) = mean(vEven(vInd_even == k)); % testing set
        mNoCV(j,k) = mean(vNoCV(vInd_noCV == k)); % no CV set
    end
    
    % find any NaNs. this is used for the normalization step because
    % mat2gray() will convert NaNs to 1. instead, we will convert NaNs to
    % 0 using these indices for visualization
    nanIndOdd  = isnan(mOdd(j,:));
    nanIndEven = isnan(mEven(j,:));
    nanIndNoCV = isnan(mNoCV(j,:));
    
    % because the sequence plot analysis is based on the MI analysis, the
    % sequence plot will be smoothed only if the MI analysis had smoothing.
    % all sequence plots are 1D so nanconv() uses 1D inputs. this smoothing
    % retains NaNs in the sequence plot
    if argins.gaussianSmoothingSigma > 0 ...
            && any(strcmp(argins.gaussianSmoothingToggle, {'both', ...
            'skaggsOnly'}))
        
        % create Gaussian filter
        imageFilter = fspecial('gaussian', [1, ...
            config.gaussianSmoothingParameters], ...
            argins.gaussianSmoothingSigma);
        
        % smooth the training, testing, and no CV sequence plots
        mOdd(j,:) = nanconv(mOdd(j,:), imageFilter, 'edge', 'nanout', '1d');
        mEven(j,:) = nanconv(mEven(j,:), imageFilter, 'edge', 'nanout', '1d');
        mNoCV(j,:) = nanconv(mNoCV(j,:), imageFilter, 'edge', 'nanout', '1d');
    end
    
    % using NaN indices from above, normalize non-NaN values in sequences
    % plots and set NaNs to zero.
    mOdd(j, ~nanIndOdd) = mat2gray(mOdd(j,~nanIndOdd)); % training set
    mOdd(j, nanIndOdd) = 0;
    
    mEven(j,~nanIndEven) = mat2gray(mEven(j,~nanIndEven)); % testing set
    mEven(j, nanIndEven) = 0;
    
    mNoCV(j, ~nanIndNoCV) = mat2gray(mNoCV(j,~nanIndNoCV)); % no CV set
    mNoCV(j,nanIndNoCV) = 0;
end
clear trialDFF tInd nanInd* i j k v*
% after ROI loop, the data for the sequence plots have been found. the next
% steps are to threshold which cells will pass the CV and sort the sequence
% plots

%% Check which ROIs pass CV and sort the sequence plots.
% remove ROIs that do not correlate well between training and testing sets
% preallocate variables
corrstore = zeros(1, length(ROIs)); % store correlation coefficients
corrstore_P = zeros(size(corrstore)); % store correlation p-values
for i = 1:length(ROIs) % loop through ROIs
    % find correlation coefficient between an ROI's values on the training
    % and testing sets
    [tempcorr, tempcorr_P] = corrcoef(mOdd(i,:), mEven(i,:));
    corrstore(i) = tempcorr(1,2); % save the correlation coefficient
    corrstore_P(i) = tempcorr_P(1,2); % save the p-value
end

% for training, testing, and no CV sets, only keep ROIs that pass alpha
% significance level specified in config
goodcorr = find(corrstore_P < config.seqPlotPvalueThreshold);
mOdd = mOdd(goodcorr,:); % training set
mEven = mEven(goodcorr,:); % testing set
mNoCV = mNoCV(goodcorr, :); % no CV set
ROIsPassedCV = ROIs(goodcorr); % save which ROIs passed the CV step

% sort the no CV sequence plot. for sorting, find which bin contains the
% max mean DFF value for each ROI. then, sort these peak values. apply
% this sorting index to the no CV matrix to sort rows by their peaks.
mNoCVunsort = mNoCV; % save unsorted sequence plot
[~, maxInd_noCV] = max(mNoCV'); %#ok<*UDIM> % find where max mean DFF is
[~, sortInd_noCV] = sort(maxInd_noCV); % sort peak DFF and save sort index
mNoCV_sort = mNoCV(sortInd_noCV, :); % apply sort index to unsorted matrix

% for the CV sequence plots, sort the training set how the non-CV sequence
% plot was sorted above. apply the sorting index from the training set to
% the testing set.
[~, maxInd_odd] = max(mOdd'); % find where max mean DFF is
[~, sortInd_odd] = sort(maxInd_odd); % sort peak DFF and save sort index

mTrain = mOdd(sortInd_odd, :); % sort training sequence plot
mTest = mEven(sortInd_odd, :); % sort testing sequence plot in the same way

% the testing sequence plot can also be sorted independently of the
% training set. by doing so, Kendall's correlation can be calculated
% between the two sorting indices derived from the training and testing
% sets individually to quantify how much they correlate.

[~, maxInd_even] = max(mEven'); % find where max mean DFF is
[~, sortInd_even] = sort(maxInd_even); % find sorting index of testing set

% calculate Kendalls correlation. check if ROIsPassedCV isn't empty. if it
% is empty, then corr() will throw an error
if ~isempty(ROIsPassedCV)
    kendallTau = corr(sortInd_odd', sortInd_even', 'type', 'kendall');
else
    kendallTau = [];
end
end