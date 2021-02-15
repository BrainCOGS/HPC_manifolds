function [output] = extractVariables(input1ROIs, input2mazeRegion, ...
    input3trials, input4activityType, input5trialLength, ...
    input6shuffle, input7eventThresholdCoefficient, ...
    input8filterdataDFF, input9imagingFilename, input10shuffleMethod, ...
    input11taskType, input12shuffleRNG, input13shuffleIndLim)
% Description: Preprocesses ROI DFF data and extracts behavioral variables
% from modeling.mat file
%
% keepTrials = find([score.trial.mainTrial]==1 & [score.trial.goodQuality]==1);
% roiGreater0o1 = find(([score.roi.transientsPerTrial]>.1)==1);
%
% Sample Call:
% out = extractVariables('all', 2, 'all', 2, 0, 0, 11, [0 0], 'E22_20170215_30p.modeling_NEW.mat', 'entireMatrix', 'towers', [], []);
% nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);
%
% ***** INPUTS *****
% input1: vector of ROI numbers to get activity
%         0 - don't extract ROI activity
%         'all' - get all ROIs activities; 'char' vector
% input2: maze regions to splice together
%         0 - cue region
%         1 - delay region
%         2 - cue + delay region
%         3 - cue + delay + trialEnd region
%         4 - entire trial (from fTrialStart to fITIEnd)
%         5 - from fTrialStart to fArmEntry
% input3: which trial numbers
%         double vector of trial numbers
%         'all' - do all trials
% input4: ROI activities to generate (digitized is boolean values for the
%         start of an ROI activation event)
%         1 - digitized event starts only
%         2 - raw activity only
%         3 - digitized activity only
% input5: double that specifies the number of observations per trial
%         upsample/downsample activity to fit a vector of given length
%
%         0 - no sampling; time data
%         !0 - sample outputs by vector length; each trial would have
%         same number of observations (resampled by position)
%
% input6: shuffle the data (circshift within trial)
%         0 - no shuffling
%         1 - shuffling
%
% input7: event threshold coefficient; input7*roi.noise = threshold for
%         event detection; applicable in binarized and event start activity types;
%         can be any number greater than zero.
%
%         5 - event threshold = 5*noise
%         10 - event threshold = 10*noise
%
% input8: filter DFF using mind_preprocess. array of two integers, first integer sets the
%       windowsize and second integer sets the thresholding
%       if [0 0] == nothing happens
%       if [0 4] == no filtering, threshold data 4*robustSTD
%       if [11 4]== filter with windowsize=11 and threshold data 4*robustSTD
%       if [11 0]== filter with windowsize=11, but don't threshold
%
% input9: character vector of the imaging session filename for an animal
%   example: 'E22_20170215_30p.modeling_NEW.mat'
%
% input10: input10shuffleMethod
%   'none' - no shuffling
%   'chunkTrials' - shuffle ROI DFF by trial. only applicable when input5 != 0
%   'entireMatrix' - circshift ROI DFF across entire table.
%   'circWithinTrials' - circshift each ROI DFF once within each trial
%
% input11: input11taskType
%   'alternation' - extract the following variables:
%       - Trial, Frame, Choice, Choice correct, Prior choice,
%       Prior choice correct, Position, Time, View angle, Yvelocity
%       Xvelocity, Velocity, Collision
%
%   'towers' - calculate the following additional variables:
%       - Evidence, EvidenceL, EvidenceR, Number of Cues, Nearest Cue Type,
%       Nearest Cue Distance
%
% input12: input12shuffleRNG
%   set RNG to this number for shuffling data
%
% input13: number of frames limit for circshifting within trials or across
% session
%
% ***** OUTPUTS *****
% output - a structure with the following fields:
%   ROIactivities - N-by-M matrix of N observations of M ROI activities
%   behavioralVariables - table of behavioral variables per observation
%   argins - structure of input arguments
%   ROIacitivitesShuffled - ROIactivities after shuffling


%% Preallocate output variables
output.ROIactivities = []; % ROI activities
output.behavioralVariables = []; % behavioral variables

%% Save inputs
argins.ROIs = input1ROIs;                   % list of ROIs
argins.mazeRegion = input2mazeRegion;       % maze region
argins.trials = input3trials;               % trial IDs
argins.activityType = input4activityType;   % ROI DFF type
argins.trialLength = input5trialLength;     % number of observations per trial
argins.eventThresholdCoefficient = input7eventThresholdCoefficient; % threshold coefficient for event detection
argins.filterdataDFF = input8filterdataDFF; % toggle to filter DFF
argins.imagingFilename = input9imagingFilename; % filename for .modeling.mat file
argins.taskType = input11taskType;          % type of task
argins.shuffle = input6shuffle;             % toggle to shuffle
argins.shuffleMethod = input10shuffleMethod; % shuffling method
argins.shuffleRNG = input12shuffleRNG;      % set shuffle RNG
argins.shuffleIndLimit = input13shuffleIndLim; % shuffle index limit

output.argins = argins;

clear input*
%% Load data (argins.imagingFilename)
D = load(argins.imagingFilename, 'score');

%% Determine which maze region (argins.mazeRegion)
% depending on the maze region, find frameStart and frameEnd for every
% trial, which will be used to extract the ROI DFF

switch argins.mazeRegion % swtich based on mazeRegion
    case 0 % cue region only
        frameStart = [D.score.trial(:).fCueEntry];
        frameEnd = [D.score.trial(:).fMemEntry];
    case 1 % delay region only
        frameStart = [D.score.trial(:).fMemEntry];
        frameEnd = [D.score.trial(:).fArmEntry];
    case 2 % cue + delay region only
        frameStart = [D.score.trial(:).fCueEntry];
        frameEnd = [D.score.trial(:).fArmEntry];
    case 3 % cue + delay + trialEnd regions
        frameStart = [D.score.trial(:).fCueEntry];
        frameEnd = [D.score.trial(:).fTrialEnd];
    case 4 % entire trial, including ITI
        frameStart = [D.score.trial(:).fTrialStart];
        frameEnd = [D.score.trial(:).fITIEnd];
    case 5 % from fTrialStart to fArmEntry
        frameStart = [D.score.trial(:).fTrialStart];
        frameEnd = [D.score.trial(:).fArmEntry];
    case 6 % manually find cue start
        for h=1:length(D.score.trial)
           tempPos = D.score.trial(h).position(:,2);
           tempiterFrame = D.score.trial(h).iterFrame;
           tempCross = find(tempPos>0);
           tempCross1 = tempCross(1);
           if abs(tempPos(tempCross1-1))<abs(tempPos(tempCross1))
               first1 = tempCross1-1;
           else
               first1 = tempCross1;
           end
           frameStart(h) = tempiterFrame(first1);
        end
        frameEnd = [D.score.trial(:).fArmEntry];
end
%% Determine which trials to use.
% argins.trials can be a character vector or a set of trial IDs (doubles).
if ischar(argins.trials) % in case that argins.trials is a character vector
    switch argins.trials % compare argins.trials to different strings
        case 'all' % do all trials
            instTrials = 1:size(D.score.trial, 2);
        case 'keepTrials' % do trials that are mainTrial and goodQuality
            instTrials = find([D.score.trial.mainTrial]==1 & ...
                [D.score.trial.goodQuality]==1);
        case 'keepTrials + leftChoice' % do keepTrials that are left choice
            instTrials = find(([D.score.trial.mainTrial]==1 & ...
                [D.score.trial.goodQuality]==1) & ...
                double([D.score.trial(:).choice]) == 1);
        case 'keepTrials + rightChoice' % do keep trials that are right choice
            instTrials = find(([D.score.trial.mainTrial]==1 & ...
                [D.score.trial.goodQuality]==1) & ...
                double([D.score.trial(:).choice]) == 2);
        case 'goodTrials' % do trials that are goodQuality
            instTrials = find([D.score.trial.goodQuality]==1);
        case 'goodTrials + leftChoice'
            instTrials = find([D.score.trial.goodQuality]==1 & ...
                double([D.score.trial(:).choice]) == 1);
        case 'goodTrials + rightChoice'
            instTrials = find([D.score.trial.goodQuality]==1 & ...
                double([D.score.trial(:).choice]) == 2);
        case 'leftChoice' % do trials that are goodQuality
            instTrials = find(double([D.score.trial(:).choice]) == 1);
        case 'rightChoice' % do trials that are goodQuality
            instTrials = find(double([D.score.trial(:).choice]) == 2);
    end
else % otherwise, argins.trials is a set of trial IDs
    instTrials = argins.trials;
end
% at this point, instTrials is a double vector of trial IDs

% using the trials specified, get the frames associated with the desired
% trials and maze region
vFrames = []; % preallocate vFrames, the vector of frame IDs
for t = instTrials % loop through trials
    % create a sequence from frameStart to frameEnd and append this
    % sequence to vFrames
    vFrames = [vFrames; [frameStart(t):frameEnd(t)]'];  %#ok<*NBRAK,*AGROW>
end

%% filter neuronal data if applicable (argins.filterdataDFF)
% filter the ROI data using mind_preprocess()
if argins.ROIs ~= 0 % only filter ROI DFF if ROIs were specified by user
    if any(argins.filterdataDFF ~= 0) % skip this if filterdataDFF is [0 0]
        dataDFF1 = D.score.dataDFF; % make copy of DFF
        % using mind_preprcoess, filter data using input arguments
        dataDFF_thresholded = mind_preprocess(dataDFF1, ...
            argins.filterdataDFF(1), argins.filterdataDFF(2));
        % replace the thresholded data into score.dataDFF
        D.score.dataDFF = dataDFF_thresholded;
    end
    
    %% extract ROI data (argins.activityType, argins.eventThresholdCoefficient)
    % extract the ROI data after the filtering step. convert the raw DFF to
    % events starts or binary towards the end of the function. this is due
    % to the event start transformation (if applicable), which must occur
    % after time or position sampling later in the function
    
    % if the ROIs specified are 'all', then all ROIs are extracted
    if ischar(argins.ROIs) && strcmp(argins.ROIs, 'all')
        argins.ROIs = 1:size(D.score.dataDFF, 2);
    end
    
    % using score.dataDFF, extract the frames and ROIs of interest. dataDFF
    % is already in terms of raw DFF activity without any event start or
    % binary transformation
    DFF = D.score.dataDFF(vFrames, argins.ROIs);
    
    % save DFF in output
    output.ROIactivities = DFF;
    
    % clear specific variables
    clear vActivity shiftInd DFFthresholded DFFdiff threshold
    clear dataDFF_filter dataDFF_thresholded dataDFF1 DFF
end

%% Extract the behavioral variables
% there are two types of behavioral variables to be extracted. the first
% are trial-wide variables, or variables that vary with each trial (Trial, MazeID,
% Choice, Choice Correct, Prior Choice, Prior Correct, Trial Probability).
% the second set of variables vary frame iterations (Position, View Angle,
% Evidence, etc.). trial-wide variables are extracted first, followed by
% frame-by-frame variables.

% preallocate the table
output.behavioralVariables = table(vFrames, 'VariableNames', {'Frame'});

% preallocate trial-wide variables
vTrials = [];               % trial ID
vMazeID = [];               % maze ID
vChoice = [];               % choice
vChoiceCorrect = [];        % correct vs error choice
vPriorChoice = [];          % prior trial's choice
vPriorCorrect = [];         % prior trial correct or error
vTrialType    = [];         % whether trial was supposed to be L or R
vTrialProb = [];            % probability of L vs R choice trial
vCollisionPos = [];            % position that mouse collide

% variables measured with iterations frequency
vPosition = [];             % y-position in maze
vPosition_X = [];           % x-position in maze
vTime = [];                 % time in seconds
vViewAngle = [];            % view angle in degrees
vYvel = [];                 % y-velocity
vXvel = [];                 % x-velocity
vVel = [];                  % velocity
vCollision = [];            % collision of mouse with walls of maze
vEvidL = [];                % leftward evidence
vEvidR = [];                % rightward evidence
vEvid = [];                 % instantaneous delta evidence (#R - #L cues)
vNumberCues = [];           % number of cues
vNCT = [];                  % nearest cue type (L or R)
vNCD = [];                  % distance between mouse and nearest cue
%vPrevEvid = [];             % evidence in previous bin/timepoint

% trial-wide variables
for t = instTrials % loop through trials
    
    % create a vector with the same number of frames for the trial
    framesVector = ones(length(frameStart(t):frameEnd(t)), 1);
    
    % for trial IDs, multiply framesVector by loop index, t
    vTrials = [vTrials; framesVector * t];
    
    % for maze IDs, multiply framesVector by the trial's maze ID
    vMazeID = [vMazeID; framesVector * D.score.trial(t).mazeID];
    
    % for choice, multiply framesVector by the animal's choice
    vChoice = [vChoice; framesVector * double(D.score.trial(t).choice)];
    
    % for choice correct, multiply framesVector by a logical vector of
    % whether the animal's choice matches the trial type
    vChoiceCorrect = [vChoiceCorrect; framesVector * (double(...
        D.score.trial(t).choice) == double(D.score.trial(t).trialType))];
    
    % for trial type, multiply framesVector by the trial's Type
    vTrialType = [vTrialType; framesVector * double(D.score.trial(t).trialType)];
    
    % for trial probability, multiply framesVector by the trial's L or R
    % probability
    vTrialProb = [vTrialProb; framesVector * double(...
        D.score.trial(t).trialProb(1))];
    
    % get the position that the mouse collided
    tempCol = find(D.score.trial(t).collision==1);
    tempPos = D.score.trial(t).position(:,2);
    tempColPos = tempPos(tempCol(1));
    vCollisionPos = [vCollisionPos; framesVector * double(tempColPos)];
    
    % for the variables describing the prior trial's choice and choice
    % correct vs error, make these NaNs for the first trial because no
    % prior trial exists
    if t == 1
        vPriorChoice = NaN(length(frameStart(t):frameEnd(t)), 1);
        vPriorCorrect = NaN(length(frameStart(t):frameEnd(t)), 1);
    else % otherwise, calculate these two variables
        vPriorChoice = [vPriorChoice; framesVector * (double(...
            D.score.trial(t-1).choice))];
        vPriorCorrect = [vPriorCorrect; framesVector * (double(...
            D.score.trial(t-1).choice) == double(...
            D.score.trial(t-1).trialType))];
    end
end % end of trial-wide variables

% variables that vary by frame iteration. several behavioral variables are
% measured in terms of 'frames' (like the ROI data). other behavioral
% variables are measured with a higher frequency with 'iterations'. to put
% the 'iteration' data on the same frequency as the 'frame' data, take the
% mean of iteration data values.
%
% some behavioral variables are only relevant for the accumulating towers
% task (including evidence). only calculate these for the towers task.
for f = 1:length(vFrames)
    trialNumber = vTrials(f); % trial ID
    % find the iteration index for this frame
    iterInd = find(vFrames(f) == D.score.trial(trialNumber).iterFrame);
    % delete indices greater than number of iterations in trial
    iterInd(iterInd > D.score.trial(trialNumber).iterations) = [];
    
    % y-position
    pos = double(mean(D.score.trial(trialNumber).position(iterInd, 2)));
    vPosition = [vPosition; pos];
    
    % x-position
    vPosition_X = [vPosition_X; double(mean(...
        D.score.trial(trialNumber).position(iterInd, 1)))];
    
    % time
    vTime = [vTime; double(mean(...
        D.score.trial(trialNumber).time(vFrames(f) == ...
        D.score.trial(trialNumber).iterFrame)))];
    
    % view angle
    vViewAngle = [vViewAngle; double(mean(...
        D.score.trial(trialNumber).position(iterInd, 3)))];
    
    % y-velocity
    yVel = double(mean(D.score.trial(trialNumber).velocity(iterInd, 2)));
    vYvel = [vYvel; yVel];
    
    % x-velocity
    xVel = double(mean(D.score.trial(trialNumber).velocity(iterInd, 1)));
    vXvel = [vXvel; xVel];
    
    % velocity
    vVel = [vVel; double(sqrt(xVel.*xVel + yVel.*yVel))];
    
    % collision
    vCollision = [vCollision; double(mean(...
        D.score.trial(trialNumber).collision(iterInd, 1)))];
    
    % variables related to the cues are only calculated in towers task
    if strcmp(argins.taskType, 'towers')
        
        % instantaneous evidence
        evidL = double(sum(vFrames(f) >= ...
            D.score.trial(trialNumber).fCueOnset{1})); % left evidence
        evidR = double(sum(vFrames(f) >= ...
            D.score.trial(trialNumber).fCueOnset{2})); % right evidence
        vEvidL = [vEvidL; evidL]; % left evidence
        vEvidR = [vEvidR; evidR]; % right evidence
        % instantaneous evidence (#R - #L cues)
        vEvid = [vEvid; evidR - evidL];
        
        % instantaneuos number of cues
        vNumberCues = [vNumberCues; evidR + evidL];
        
        % calculate neareast cue distance / nearest cue type
        % distance to nearest L cue
        distL = double(min(abs(pos - ...
            D.score.trial(trialNumber).cuePos{1})));
        % distance to nearest R cue
        distR = double(min(abs(pos - ...
            D.score.trial(trialNumber).cuePos{2})));
        
        % the closest cue is left because either distR is empty (for trials
        % with only left cues, for example), or distL is less than distR.
        if isempty(distR) % if no rightward cue
            if isnan(distL) % and no leftward cue, append NaNs
                vNCT = [vNCT; nan];
                vNCD = [vNCD; nan];
            else % if leftward cue exists, then append leftward values
                vNCT = [vNCT; 1];
                vNCD = [vNCD; distL];
            end
        end
        if distL < distR % distance to L cue is less than distance to R cue
            vNCT = [vNCT; 1];
            vNCD = [vNCD; distL];
        end
        
        % otherwise, the closest cue is right because either distL is empty
        % (for trials with only R cues), or distR is less than distL.
        if isempty(distL) % no L cue
            if isnan(distR) % no R cue
                vNCT = [vNCT; nan]; % append NaNs
                vNCD = [vNCD; nan];
            else % otherwise, append R values
                vNCT = [vNCT; 2];
                vNCD = [vNCD; distR];
            end
        end
        if distR < distL % distance to R cue is less than distance to L cue
            vNCT = [vNCT; 2];
            vNCD = [vNCD; distR];
        end
        
        % if both distL and distR are NaNs, then append NaNs
        if isnan(distL) & isnan(distR)  %#ok<*AND2>
            vNCT = [vNCT; nan];
            vNCD = [vNCD; nan];
        end
    end
end
clear distL distR evidL evidR frame* iterInd pos t trialNumber f xVel yVel pos_X

% subtract 1 from choice and NCT to make binary (0 or 1), not (1 or 2)
vChoice = vChoice - 1;
vTrialType = vTrialType -1;
vPriorChoice = vPriorChoice - 1;
vNCT = vNCT - 1;

% save behavioral variables
output.behavioralVariables.Trial = vTrials;
output.behavioralVariables.MazeID = vMazeID;
output.behavioralVariables.Choice = vChoice;
output.behavioralVariables.ChoiceCorrect = vChoiceCorrect;
output.behavioralVariables.TrialType = vTrialType;
output.behavioralVariables.PriorChoice = vPriorChoice;
output.behavioralVariables.PriorCorrect = vPriorCorrect;
output.behavioralVariables.TrialProb = vTrialProb;
output.behavioralVariables.Position = vPosition;
output.behavioralVariables.Position_X = vPosition_X;
output.behavioralVariables.Time = vTime;
output.behavioralVariables.ViewAngle = vViewAngle;
output.behavioralVariables.Yvelocity = vYvel;
output.behavioralVariables.Xvelocity = vXvel;
output.behavioralVariables.Velocity = vVel;
output.behavioralVariables.Collision = vCollision;
output.behavioralVariables.CollisionPos = vCollisionPos;

% for towers task, save additional variables
if strcmp(argins.taskType, 'towers')
    output.behavioralVariables.EvidenceL = vEvidL;
    output.behavioralVariables.EvidenceR = vEvidR;
    output.behavioralVariables.Evidence = vEvid;
    output.behavioralVariables.EvidenceSmooth = smooth(vEvid, 11);
    output.behavioralVariables.NumberOfCues = vNumberCues;
    output.behavioralVariables.NearestCueType = vNCT;
    output.behavioralVariables.NearestCueDistance = vNCD;
    
    % Add in the previous evidence, which is -1 or 1 or NaN
    %vEvidPrev = nan(length(vEvid),1);
    vEvidPrev = [];
    
    % Need to split into trials because of the trial starts
    trialn = unique(vTrials);
    
    for j=1:length(trialn)
        vEvid_trial = vEvid(vTrials==trialn(j));
        vtrialType_trial = vTrialType(vTrials==trialn(j));
        vEvidAbs = abs(vEvid_trial);
        vEvid1 = vEvidAbs;
        vEvid1 = [0; vEvid1];
        vEvid1_diff = diff([vEvid1]);
        vEvid1_ind = find(vEvid1_diff~=0);
        vEvidPrev_trial = nan(length(vEvid_trial),1);
        for i=1:length(vEvid1_ind)
            if i~=length(vEvid1_ind)
                if vEvid1(vEvid1_ind(i)+1)>vEvid1(vEvid1_ind(i))
                    vEvidPrev_trial(vEvid1_ind(i)+1:vEvid1_ind(i+1)) = 1;
                else
                    vEvidPrev_trial(vEvid1_ind(i)+1:vEvid1_ind(i+1)) = -1;
                   % vEvid1 = -1*vEvid1;
                end
            else
                if vEvid1(vEvid1_ind(i)+1)>vEvid1(vEvid1_ind(i))
                    vEvidPrev_trial(vEvid1_ind(i)+1:length(vEvid1)) = 1;
                else
                    vEvidPrev_trial(vEvid1_ind(i)+1:length(vEvid1)) = -1;
                   % vEvid1 = -1*vEvid1;
                end
            end
        end
        % have to remove the effect of that zero
        vEvidPrev_trial = vEvidPrev_trial(2:end);
        vEvidPrev = [vEvidPrev; vEvidPrev_trial];
    end
    
    output.behavioralVariables.PrevEvidence = vEvidPrev;
    
end
clear v*

%% sample each trial to specified length (input5activityVectorLength)
% if argins.trialLength is not zero, then convert the ROI and DFF data from
% sampling by position instead of time. if it is zero, then retain time
% sampling (by frames).
if argins.trialLength ~= 0 % if 0, don't sample
    % determine task type because towers and alternation tasks have
    % different maze region lengths
    switch argins.taskType
        case {'towers','alternation'}
            switch argins.mazeRegion
                case 0 % cue region
                    posStart = 0;
                    posEnd = 200;
                case 1 % delay region
                    posStart = 200;
                    posEnd = 300;
                case 2 % cue + delay region
                    posStart = 0;
                    posEnd = 300;
                case 3 % cue + delay + trialEnd
                    posStart = 0;
                    posEnd = 300;
                case 4 % fTrialStart to fITIEnd
                    posStart = -30;
                    posEnd = 300;
                case 5  % from fTrialStart to fArmEntry
                    posStart = -30;
                    posEnd = 300;
            end
        case {'AlternationJeff','alternationJeff'}
            switch argins.mazeRegion
                case 4 % fTrialStart to fITIEnd
                    posStart = -30;
                    posEnd = 348;
                case 6
                    posStart = -30;
                    posEnd = 348;
            end
    end
    % using start and end positions, create a uniform position sequence
    ratio = (posEnd - posStart)/argins.trialLength;
    seq = posStart:ratio:posEnd;
    
    % get position and trial IDs
    vPos = output.behavioralVariables.Position;
    vTrials = output.behavioralVariables.Trial;
    
    % for each value of seq, find the closest corresponding frame in the
    % time sampling. ind will be indexed into the data sampled by time
    % (frames).
    ind = zeros(1, length(seq)*length(instTrials));
    counter = 1;
    for t = instTrials % loop through trials
        ind2 = find(vTrials == t); % get relevant trial observations
        vPos2 = NaN(size(vPos)); % preallocate a position vecotr
        vPos2(ind2) = vPos(ind2); % copy only pos values for this trial
        for i = 1:length(seq) % going through seq positions
            % find the vPos2 position that is closest to the sequence value
            [~, ind(counter)] = min(abs(vPos2 - seq(i)));
            counter = counter + 1;
        end
    end
    
    % use ind on the original data. depending on the value of
    % input5activityVectorLength, this may delete frames of data (if the
    % input is less than the number of frames in a trial), or it may copy
    % frames of data (if the input is greater than the number of frames in
    % a trial). after this step, all trials have the same number of
    % observations.
    output.ROIactivities = output.ROIactivities(ind, :);
    output.behavioralVariables = output.behavioralVariables(ind, :);
    output.positionBins = seq;
end

clear counter i ind* v* ratio seq t
clear xVel yVel

%% Trasnform ROI DFF activity (if applicable) (input4activityType)
% there are two types of transformations for the DFF data: binarization and
% event starts. binarization sets values above a threshold to 1 and values
% below to 0. event starts also binarizes the data but only allows the
% first value of an activation event that passes the threshold to be one
% (the start of the event). all other values are zero.

DFF = output.ROIactivities; % copy DFF data

noiseROI = zeros(1,size(DFF,2));
for i=1:size(DFF,2)
    noiseROI(i) = robustSTD(DFF(:,i));
end

switch argins.activityType
    case 1 % digitized event starts only
        % threshold DFF data using input argument and the ROI's noise value
        DFFthresholded = DFF > (argins.eventThresholdCoefficient * ...
            noiseROI);
        %    [D.score.roi(argins.ROIs).noise]);
        % take diff to find event starts where diff is 1
        DFFdiff = diff(DFFthresholded);
        DFFdiff(DFFdiff == -1) = 0;
        % diff loses a row, so use the first row of DFFthresholded as the
        % first row.
        DFF = [DFFthresholded(1,:); DFFdiff];
    case 3 % digitzed activity only
        % threshold DFF data, similar to event starts
        DFFthresholded = DFF > (argins.eventThresholdCoefficient*...
            noiseROI);
            %[D.score.roi(argins.ROIs).noise]);
        DFF = DFFthresholded;
end

% put DFF back into the output
output.ROIactivities = DFF;

%% Shuffle neuronal activities if applicable (argins.shuffle)
% shuffle the data using the specified shuffling method
if argins.shuffle == 1
    rng(argins.shuffleRNG); % set shuffle RNG
    % preallocate shuffled DFF activity
    output.ROIactivitiesShuffled = zeros(size(output.ROIactivities));
    switch argins.shuffleMethod
        
        case 'none'
            disp('ERROR: NEED A SHUFFLE METHOD');
        
        % circular shift within each trial
        case 'circWithinTrials'
            output.ROIactivitiesShuffled = circshift_within_trials(...
                output.ROIactivities, output.behavioralVariables.Trial, ...
                argins.shuffleIndLimit);
            
            % trial-by-trial shuffling per ROI; preserves ROI activity
            % within a trial, but shuffles chunks of trials together
        case 'chunkTrials'
            % this shuffling method only works if each trial has the same
            % number of observations (argins.trialLength ~= 0). if trials
            % are not the same length, then shuffling chunks of them
            % independently for each ROI will not preserve trial identity
            % across ROIs.
            if argins.trialLength ~= 0
                % preallocate trial indices
                output.ROIactivitiesShuffledTrialIndices = ...
                    zeros(length(unique(...
                    output.behavioralVariables.Trial)), ...
                    size(output.ROIactivities, 2));
                
                for i = 1:size(output.ROIactivities, 2) % loop through ROIs
                    vals = output.ROIactivities(:,i); % get DFF values
                    
                    % get trials and trialIDs
                    vTrials = output.behavioralVariables.Trial;
                    t = unique(vTrials);
                    
                    % shuffle trialIDs
                    tShuffle = t(randperm(length(t)));
                    % preallocate shuffle vals
                    valsShuffle = vals;
                    for ii = 1:length(t) % loop through trials
                        % assign real vals with shuffle trial IDs to
                        % shuffle value set
                        valsShuffle(vTrials == t(ii)) = ...
                            vals(vTrials == tShuffle(ii));
                    end
                    % save the results
                    output.ROIactivitiesShuffled(:, i) = valsShuffle;
                    output.ROIactivitiesShuffledTrialIndices(:, i) = ...
                        tShuffle;
                end
            else
                error(['Trial lengths must be the same for this ' ...
                    'shuffle method.']);
            end
            
            % circshift each ROI's DFF throughout entire DFF matrix.
            % this does not take the trials into account
        case 'entireMatrix'
            output.ROIactivitiesShuffled = ...
                zeros(size(output.ROIactivities));
            for i = 1:size(output.ROIactivities, 2) % loop through ROIs
                vals = output.ROIactivities(:,i); % get DFF values
                vals = circshift(vals, ... % do a circshift
                    randi([1+argins.shuffleIndLimit, ...
                    length(vals) - argins.shuffleIndLimit]), 1);
                output.ROIactivitiesShuffled(:, i) = vals; % save result
            end
    end
end
end