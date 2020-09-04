%% Function to generate the cross-validation set
function CV = generateCrossValSet_v2(behavioralVariables, numFolds)
% Description: create a cross-validation set based on the table in
% behavioralVariables. because temporal correlations exist in both animal
% behavior and ROI activity for a frame and its local neighbors, a random
% assignment procedure for partitioning the data into cross-validation sets
% is not appropriate. this random assignment procedure may allow frame (t) 
% and frame (t+1) to be split into different training sets, which may yield
% an overly optimistic testing accuracy and be less conservative. instead,
% this cross-validation procedure randomly assigns entire trials, in which
% local temporal correlations are contained, to different training and
% testing sets. 
%
% Sample call:
% CV = generateCrossValSet_v2(behavioralVariables, numFolds);
%
% ***** INPUTS *****
% behavioralVariables - table out from extractVariables()
% numFolds - number of k-folds to split the data into
% 
% ***** OUTPUTS *****
% CV - a 1-by-K structure with multiple fields, where K is number of folds:
%       test            - trial IDs in test set
%       train           - trial IDs in train set
%       testLocations   - rows of m in test set
%       trainLocations  - rows of m in train set
%       testData        - m with testLocations only
%       trainData       - m with trainLocations only
%       specs           - function specifications

% preallocate output structure
CV(numFolds).test = [];

trials      = behavioralVariables.Trial; % trial IDs for each observation
uniqTrials  = unique(trials);            % list of trial IDs
numTrials   = length(uniqTrials);        % number of trials
splitTrials = floor(numTrials/numFolds); % number of trials per split
vFolds      = repelem(1:numFolds, splitTrials)'; % fold index per trial

% because splitTrials may round, vFolds may not have a fold index for every
% trial. append more values to vFolds in this case
if size(vFolds, 1) ~= numTrials 
    sizeDiff = abs(size(vFolds, 1) - numTrials); % size difference 
    vFolds = [vFolds; (1:sizeDiff)']; % append this sequence to vFolds
end
clear sizeDiff

% randomize vFolds
vFolds = vFolds(randperm(length(vFolds)));

% break the data into different folds
for i = 1:numFolds % loop through folds
    CV(i).test              = uniqTrials(vFolds == i); % test trial IDs
    CV(i).train             = uniqTrials(vFolds ~= i); % train trial IDs
    
    % observations in test set
    CV(i).testLocations     = find(ismember(trials,CV(i).test)==1);
    % observations in train set
    CV(i).trainLocations    = find(ismember(trials,CV(i).train)==1);
    % test data
    CV(i).testData          = behavioralVariables(CV(i).testLocations,:);
    % train data
    CV(i).trainData         = behavioralVariables(CV(i).trainLocations,:);
end
% save input arguments
CV(1).specs.behavorialVariables = behavioralVariables;
CV(1).specs.numFolds = numFolds;
end