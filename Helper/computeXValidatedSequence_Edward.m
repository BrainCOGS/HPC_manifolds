% COMPUTEXVALIDATEDSEQUENCE       Definition of choice-specific sequences a la Harvey 2012, with
%                                 adjustments for GCamp6f data and additional statistical tests.
%
%  Inputs
% ========
%   modelFile       : .modeling.mat file as output from computeActivityResponse().
%   lazy            : If true, doesn't do anything if the output file already exists.
%   versionInfo     : Optional precomputed code versioning info [output from collectVersionInfo()]
%                     to speed up execution when calling this function sequentially on many
%                     modelFile inputs. If an empty vector [] is provided, will call
%                     collectVersionInfo() on its own.
%   trialSetName    : What set of trials (stored in the modelFile) to use for this computation. If
%                     not provided or is an empty vector, this code will use 'trialSet' which
%                     corresponds to main maze trials. Can be set to 'guidedSet' instead to use
%                     visually guided maze trials.
%
%  Output
% ========
%   Creates a .sequences.mat file in the same directory as the input modelFile, with the following
%   contents:
%     versionInfo               : code version and data sources; the analysis configuration cfg in 
%                                 this code is stored as analysisParams in this output struct
%     globalID                  : IDs of all selected cell (with enough transients/trial; see cfg)
%     trialSet                  : indices of trialType-by-outcome trials in epochDFF below. Rows of
%                                 this matrix correspond to left and right rewarded trials in that
%                                 order. Columns of this matrix correspond to correct and error
%                                 outcomes of the animal's choice. i.e. Access as:
%                                     trialSet{L/R rewarded,correct/error}
%     
%     sortTrials                : Trials used for sorting: correct, odd-indexed trials only, i.e.:
%                                     sortTrials{L/R/all choices}
%     testTrials                : Trials left out of the sortTrials set, i.e. correct trials if
%                                 even-indexed (first column), or any error trials (second column)
%     trialIDs                  : original IDs of trials in trialSet.
%     epochDuration             : duration of each task epoch
%     epochDFF                  : time-by-trials-by-cell matrix of dF/F of cells. Here time is
%                                 defined as discrete epoch bins so that all trials can be
%                                 registered onto the same time axis. This is defined via a
%                                 piecewise linear mapping between actual wall clock time during the
%                                 trial and an integer numbering of the task epochs (start, cue
%                                 region, delay period, turn region, ITI).
%     selectivity               : 
%     preferredSide             : 
%     epochCenters              : centers of epoch bins used to compute the bin-averaged dF/F
%     peakFrame                 : 
%     peakEpoch                 : 
%     firingFieldRange          : 
%     cellOrder                 : 
%     ridgeToBkg                : 
%     reliability               : 
%     reliabilityCI             : 
%     epochBinNoise             : 
%     shuffledRidge             : 
%     ridgeSignificance         : significance of the ridge-to-background ratio, where the null
%                                 hypothesis is formed by randomly rotating the data in time.
%   
function temp3a = computeXValidatedSequence(modelFile, input7eventThresholdCoefficient, input8filterdataDFF, numWorkers, lazy, versionInfo, trialSetName)
  
argins.eventThresholdCoefficient = input7eventThresholdCoefficient;
argins.filterdataDFF             = input8filterdataDFF;

  %% Default arguments
  if exist(modelFile,'file') ~= 2
    candidateInput        = rdir(fullfile(modelFile, '*.modeling.mat'));
    modelFile             = candidateInput.name;
  end
  
  if nargin > 3 && ~isempty(numWorkers)
    startParallelPool(-numWorkers);
  end
  
  if nargin < 5 || isempty(lazy)
    lazy                  = false;
  end
  if nargin < 6
    versionInfo           = [];
  end
  if nargin < 7 || isempty(trialSetName)
    trialSetName          = 'trialSet';
  end
  
  %% Analysis configuration
  [inputDir,inputName]    = parsePath(modelFile);
  cfg.name                = strrep(regexprep(inputName, '(um)_.*', '$1'), '_', ' ');
  cfg.outputFile          = fullfile(inputDir, regexprep(inputName, '[.].*', '.sequences-xval.mat'));

  cfg.minTransPerTrial    = 0.1;
  cfg.minActiveFrac       = 0.25;
  cfg.minActiveFrames     = 2;
  cfg.noiseFactor         = 3:5;
  cfg.minFieldHeight      = 0.5;

  cfg.behaviorVar         = {'position', 'velocity', 'sensorVel', 'sumCues', 'leakyNCues', 'basisCues', 'basisCues2', 'cueAtLag'};
  cfg.rewardEpoch         = 4;
%   cfg.integralTau         = 0.2:0.2:4;   % seconds
  cfg.integralTau         = [0.2:0.2:1, 1.25:0.25:2, 2.5:0.5:3 4:5];   % seconds
  cfg.integralResolution  = 1e-3;
  cfg.basisScale          = 1/4;
  cfg.basisTau            = exp(linspace(log(0.1), log(6), 10));
%   cfg.basisTau2           = linspace(0.1^(1/2), 6^(1/2), 15).^2;
  cfg.basisTau2           = linspace(0.1^(1/2), 6^(1/2), 10).^2;
  cfg.cueWindow           = 0.3;                        % seconds
  cfg.cuePatternBins      = 2:5;
  
%   cfg.epochBinning        = [8 17 12 5 11 28];          % roughly 200ms
  cfg.epochBinning        = [0 41 21 0 0 0];          % roughly 250ms, from mean epochDuration across all data
  cfg.epochEdges          = accumfun(2, @(x,y,z) butlast(linspace(x,y,z)), 0:numel(cfg.epochBinning)-1, 1:numel(cfg.epochBinning), cfg.epochBinning);
  cfg.epochSel            = cfg.epochEdges < 5;
  cfg.numShuffles         = 1000;
  cfg.numSides            = numel(Choice.all());


  %% Load data
  if lazy && exist(cfg.outputFile, 'file')
    fprintf('----  %s\n', modelFile);
    return;
  else
    fprintf('****  %s\n', modelFile);
  end
  
  load(modelFile, '-regexp', '^(?!cfg$).*');

  noData                  = arrayfun(@(x) isempty(x.trial), score);
  score(noData)           = [];
  event(noData)           = [];
  assert( numel(score) == 1 );
  
  cfg.deltaT              = score.deltaT;
  
  %% Version info tracking
  if ~exist('collectVersionInfo','file') || PaneledFigure.developmentMode() || exist('output', 'var')
    versionInfo           = struct('parameters', cfg);
  elseif isempty(versionInfo)
    versionInfo           = collectVersionInfo([], cfg, [], [], modelFile);
  else
    versionInfo           = collectVersionInfo(versionInfo, cfg, [], [], modelFile);
  end
  
  
  %%  || Additional processing of behavioral data
  %   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  %% Define leaky (exponential) integration kernels
  tauFrames               = cfg.integralTau / cfg.deltaT;
  maxFrames               = ceil(expinv(1-cfg.integralResolution, tauFrames));
  intgKernel              = arrayfun(@(x,y) differenceExponentials(0.5:x,0,y), maxFrames, tauFrames, 'UniformOutput', false);
  cueWinKernel            = ones(round( cfg.cueWindow / cfg.deltaT ),1);
  tauLagIndex             = round(tauFrames);
  
  %% Define bump-basis kernels
  tauFrames               = cfg.basisTau / cfg.deltaT;
  basisWidth              = tauFrames * cfg.basisScale;
  basisK                  = ( tauFrames ./ basisWidth ).^2;
  basisTheta              = basisWidth.^2 ./ tauFrames;
  maxFrames               = ceil(gaminv(1-cfg.integralResolution, basisK, basisTheta));
  basisKernel             = arrayfun(@(x,y,z) discretizedGamma(0:x,y,z), maxFrames, basisK, basisTheta, 'UniformOutput', false);

  %%
  tauFrames               = cfg.basisTau2 / cfg.deltaT;
  basisWidth              = tauFrames.^(0.5);
  basisK                  = ( tauFrames ./ basisWidth ).^2;
  basisTheta              = basisWidth.^2 ./ tauFrames;
  maxFrames               = ceil(gaminv(1-cfg.integralResolution, basisK, basisTheta));
  basisKernel2            = arrayfun(@(x,y,z) discretizedGamma(0:x,y,z), maxFrames, basisK, basisTheta, 'UniformOutput', false);
%   longfigure; plot((0:max(maxFrames)) * cfg.deltaT, catpad(2, 0, basisKernel2{:}))
  
  %% Per-trial processing
  event.leakyNCues        = zeros([size(event.sumCues), numel(tauLagIndex)]);
  event.basisCues         = zeros([size(event.sumCues), numel(basisKernel)]);
  event.basisCues2        = zeros([size(event.sumCues), numel(basisKernel2)]);
  event.cueAtLag          = false([size(event.sumCues), numel(tauLagIndex)]);
  for iTrial = 1:numel(score.trial)
    trial                 = score.trial(iTrial);
    itiFrames             = trial.fMemEntry:trial.fITIEnd;
    allFrames             = trial.fCueEntry:trial.fITIEnd;
    for iSide = 1:cfg.numSides
      %% Define cue information to be constant past the end of the cue period
      event.sumCues(itiFrames,iSide) = trial.numCues(iSide);
      if trial.numCues(iSide) < 1
        continue;
      end
      
      %% Compute various filtered versions of cue onsets
      cueOnset            = event.cueOnset(allFrames,iSide);
      cueInWindow         = filter(cueWinKernel, 1, cueOnset(:), [], 1) > 0;
      cueInWindow         = padarray(cueInWindow, [tauLagIndex(end),0], false, 'post');
      for iTau = 1:numel(tauLagIndex)
        %% Convolutions of pulses with exponential kernels
        convCues          = conv(cueOnset, intgKernel{iTau});
        event.leakyNCues(allFrames,iSide,iTau)  = convCues(1:numel(allFrames));
        
        %% Is there a cue within a particular time window, that occurred at a given time lag ago?
        event.cueAtLag(allFrames(tauLagIndex(iTau):end),iSide,iTau)   ...
                          = cueInWindow(1:numel(allFrames) - tauLagIndex(iTau) + 1);
      end
      
      %% Convolutions of pulses with bump-basis kernels
      for iTau = 1:numel(basisKernel)
        convCues          = conv(cueOnset, basisKernel{iTau});
        event.basisCues(allFrames,iSide,iTau)  = convCues(1:numel(allFrames));
      end
      
      %% Convolutions of pulses with bump-basis kernels
      for iTau = 1:numel(basisKernel2)
        convCues          = conv(cueOnset, basisKernel2{iTau});
        event.basisCues2(allFrames,iSide,iTau) = convCues(1:numel(allFrames));
      end
    end
  end
  
  %% Convert raw Arduino sensor velocity (dots/CPU-tick) to cm/s -- N.B. this is approximate because the actual conversion factor, although stored, is not conveniently available
  % Also two things can be wrong about this calculation:
  %   (1) In some versions of the protocol the view angle is not restricted to 0 in the start period
  %   (2) Some movement rules have a nonlinear relationship between sensor velocity and
  %       virtual-world velocity
  % However for now I'm fine with these issues since the only intention is to convert sensor
  % dots/sec to a physical scale cm/s, which does not have to be perfect.
  
  selPreCue               = event.epoch < 0.5;     % here the mouse is restricted to travel straight down the maze, so virtual y velocity and forward velocity of the ball are exactly equivalent
  if std(event.velocity(selPreCue,1)) > 1e-5
    warning ( 'computeXValidatedSequence:sensorVel', 'Nonzero virtual x-velocity in pre-cue period (std. dev. = %.3g).'   ...
            , std(event.velocity(selPreCue,1)) );
  end
  
  virtualVel              = event.velocity(selPreCue,2);
  sensorVel               = event.sensorVel(selPreCue,2);
  sel                     = ~isnan(virtualVel) & ~isnan(sensorVel);
  velocityScale           = sensorVel(sel) \ virtualVel(sel);
  if ~isfinite(velocityScale)
    error ( 'computeXValidatedSequence:sensorVel'                                                           ...
          , 'Invalid scale factor %.3g/%.3g computed for converting raw Arduino sensor velocity to cm/s.'   ...
          , mean(virtualVel(sel), 'omitnan'), mean(sensorVel(sel), 'omitnan') );
  end
	event.sensorVel(:,1:2)  = event.sensorVel(:,1:2) * velocityScale;
      
  
  %%  ||   Basic cell/trial selection, and time warping to behavioral epochs
  %   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  %% Select cells with a high enough activity rate for analysis
  activeCells             = [score.roi.transientsPerTrial] >= cfg.minTransPerTrial;
  globalID                = [score.roi(activeCells).globalID];

  %% Make sure that all selected trials have mazes with the same cue region length
  trialIDs                = score.(trialSetName);           % ordered set of { correct L-rewarded , error L-rewarded
  selTrials               = score.trial([trialIDs{:}]);     %                ; correct R-rewarded , error R-rewarded } trials
  if ~isempty(selTrials)
    mazeConfig            = [selTrials.mazeConfig];
    mazeLength            = [mazeConfig.lCue; mazeConfig.lMemory]';
    [mazeLength,~,iMaze]  = unique(mazeLength, 'rows');
    if size(mazeLength,1) > 1
      %% Keep only trials with the most common maze type
      nMazeTrials         = aggregateByBin(ones(size(iMaze)), iMaze, size(mazeLength,1), 1, @sum);
      [~,selMaze]         = max(nMazeTrials);
      sel                 = iMaze == selMaze;
      warning ( 'computeXValidatedSequence:mazeLength'                                                          ...
              , '%d different types of mazes found in %s. Keeping only %d/%d trials belonging to maze type %d.'   ...
              , size(mazeLength,1), trialSetName, sum(sel), numel(iMaze), selMaze                                 ...
              );

      selTrials           = cellfun(@(x,y) y(sel(x)), span2index(cellfun(@numel, trialIDs(:))), trialIDs(:), 'UniformOutput', false);
      trialIDs            = reshape(selTrials, size(trialIDs));
    end
  end
  
  %% Select subset of data corresponding to main maze trials; note that this can generate a small number of discontinuities 
  selTrials               = score.trial(sort([trialIDs{:}]));     % ordered set as in trialIDs
  if isempty(selTrials)
    warning('computeXValidatedSequence:selTrials', 'Dataset does not contain any valid trials in "%s".', trialSetName);
    return;
  end
  
  selFrames               = [selTrials.frame];
  cellDFF                 = score.dataDFF(selFrames,:);
  cellDFF                 = mind_preprocess(cellDFF, argins.filterdataDFF(1), argins.filterdataDFF(2));
  epoch                   = event.epoch(selFrames,:);
  trialFrames             = span2index([selTrials.numCNMF]);
  
  % Subsets of trials: trialSet{L/R rewarded,correct/error}
  trialSet                = span2index(cellfun(@numel, trialIDs(:)));
  trialSet                = reshape(trialSet, size(trialIDs));
  
  % Trials used for sorting: correct, odd-indexed trials -> sortTrials{L/R/all choices}
  sortTrials              = cellfun(@(x) x(mod(x,2)==1), trialSet(:,1), 'UniformOutput', false);
  sortTrials{end+1}       = [sortTrials{1:cfg.numSides}];
  testTrials              = cellfun(@(x) x(mod(x,2)==0), trialSet(:,1), 'UniformOutput', false);
  testTrials(:,end+1)     = trialSet(:,end);
  testTrials(end+1,:)     = catcell(2, testTrials, 1);
  
  %% Trial-based behavioral variables
%   behavior                = struct();
%   behavior.cuePattern     = cell(size(cfg.cuePatternBins));
%   for iPattern = 1:numel(cfg.cuePatternBins)
%     %% Count number of cues within a given division of the cue-region epoch
%     patternEpoch          = linspace(1, 2, cfg.cuePatternBins(iPattern) + 1);
%     behavior.cuePattern{iPattern} = nan(numel(selTrials), cfg.numSides, cfg.cuePatternBins(iPattern));
%     
%     for iTrial = 1:numel(selTrials)
%       trial               = selTrials(iTrial);
%       for iSide = 1:cfg.numSides
%         cueEpoch          = event.epoch(trial.fCueOnset{iSide});
%         behavior.cuePattern{iPattern}(iTrial,iSide,:)       ...
%                           = histfast(cueEpoch, patternEpoch, [], true, false);
%       end
%     end
%   end
  
  %% Keep track of repeated trials
%   [repeatedTrials, trialMatchInfo]    ...
%                           = groupRepeatedTrials(selTrials);
%   behavior.repeatedTrials = cellfun(@uint16, repeatedTrials, 'UniformOutput', false);

  %% Record typical duration of each epoch and decide parameters for averaging data in epoch bins 
  totEpochs               = round(max(event.epoch(:)));
  refEpoch                = 0:totEpochs;
  epochDur                = nan(numel(trialFrames), totEpochs);
  for iTrial = 1:numel(trialFrames)
    assert(issorted(epoch(trialFrames{iTrial})));
    
    iEpoch                = binarySearch(double(epoch(trialFrames{iTrial})), refEpoch, 0, 2);
    epochDur(iTrial,:)    = diff(iEpoch);
  end
  epochDur(epochDur == 0) = nan;
  epochDuration           = mean(epochDur, 1, 'omitnan') * cfg.deltaT;
  epochBin                = binarySearch(cfg.epochEdges, epoch, -1, -1);
  
  %% Collect (epoch,trial,cell) dF/F for sequence computations 
  % N.B. There can be NaNs in epochDFF and associated epoch-based behavioral quantities, because in
  % some trials two imaging frames could correspond ambiguously to one epoch bin 
  
  dataTrial               = accumfun(1, @(x,y) repmat(x,y,1), 1:numel(selTrials), [selTrials.numCNMF]);
  dataBins                = [epochBin(:), dataTrial(:)];
  nDataBins               = [numel(cfg.epochEdges), numel(selTrials)];
  epochDFF                = averageInBin(cellDFF, dataBins, nDataBins, 1);
  
  for i=1:size(epochDFF,3)
      for j=1:size(epochDFF,2)
          
          tempDFF = epochDFF(:,j,i);
          while sum(isnan(tempDFF))~=0
              nanDFF = find(isnan(tempDFF));
              tempDFF(nanDFF) = tempDFF(nanDFF-1);
          end
          epochDFF1(:,j,i) = tempDFF;          
      end
  end
  
  for i=1:size(epochDFF1,3)
     DFF1 = epochDFF1(:,:,i); 
     DFF(:,i) = DFF1(:);
  end
  
  % Copied and pasted this section from extractVariables
noiseROI = zeros(1,size(DFF,2));
for i=1:size(DFF,2)
    noiseROI(i) = robustSTD(DFF(:,i));
end

DFFthresholded = DFF > (argins.eventThresholdCoefficient * ...
    noiseROI);
%    [D.score.roi(argins.ROIs).noise]);
% take diff to find event starts where diff is 1
DFFdiff = diff(DFFthresholded);
DFFdiff(DFFdiff == -1) = 0;
% diff loses a row, so use the first row of DFFthresholded as the
% first row.
DFF = [DFFthresholded(1,:); DFFdiff];

epochDFF2 = reshape(DFF, 60, size(epochDFF,2), []);

  for i=1:size(epochDFF,2)
     temp3a{i} = squeeze(epochDFF2(:,i,:)); 
  end
  
  %% Store some behavioral quantities, similarly averaged in the epoch bins
%   for field = cfg.behaviorVar
%     variable              = event.(field{:})(selFrames,:);
%     behavior.(field{:})   = averageInBin(variable, dataBins, nDataBins, 1);
%     if islogical(variable)
%       behavior.(field{:}) = behavior.(field{:}) > 0;
%     else
%       behavior.(field{:}) = single(behavior.(field{:}));
%     end
%     
%     dataSize              = size(event.(field{:}));
%     if numel(dataSize) > 2
%       behavior.(field{:}) = reshape(behavior.(field{:}), [size(behavior.(field{:}),1), size(behavior.(field{:}),2), dataSize(2:end)]);
%     end
%   end
%   
%   %% Set position variables to constant during the ITI
%   iTrialEnd               = find(cfg.epochEdges >= cfg.rewardEpoch, 1, 'first');
%   behavior.position(iTrialEnd:end,:,:)   ...
%                           = repmat(shiftdim(accumfun(1, @(x) x.position(end,:), selTrials), -1), numel(cfg.epochEdges)-iTrialEnd+1, 1, 1);
% 
%   
%   %% Past-trial behavioral info, or NaN for the first trial in the session
%   allTrialIDs             = [trialIDs{:}];
%   hasPast                 = allTrialIDs > 1;
%   pastTrials              = score.trial(allTrialIDs(hasPast) - 1);
%   
%   behavior.pastChoice     = nan(1,numel(selTrials), 'single');
%   behavior.pastReward     = nan(1,numel(selTrials), 'single');
%   behavior.pastSumCues    = nan(1,numel(selTrials),size(behavior.sumCues,3), 'single');
%   behavior.pastPosition   = nan(1,numel(selTrials),size(behavior.position,3), 'single');
%   behavior.pastChoice(hasPast)        = [pastTrials.choice];
%   behavior.pastReward(hasPast)        = [pastTrials.choice] == [pastTrials.trialType];
%   behavior.pastSumCues(:,hasPast,:)   = shiftdim(cat(1, pastTrials.numCues), -1);
%   behavior.pastPosition(:,hasPast,:)  = shiftdim(accumfun(1, @(x) x.position(end,:), pastTrials), -1);
%   
%   %% HACK for adding to behavioral data storage
%   %{
%   if exist('output', 'var')
%     output.versionInfo.parameters     = cfg;
%     for what = fieldnames(behavior)'
%       if isfield(output, what{:})
%         assert(isequaln( behavior.(what{:}), output.(what{:}) ));
%       else
%         output.(what{:})  = behavior.(what{:});
%       end
%     end
%     
%     save( cfg.outputFile, '-struct', 'output');
%     return;
%   end
%   %}
%   
%   
%   %%  ||   Selection criteria for consistently active cells with localized activity
%   %   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%   %% For each cell, estimate the noise for epoch-binned dF/F by assuming a median number of frames per epoch bin
%   framesPerEpochBin       = bsxfun(@rdivide, epochDur, cfg.epochBinning);
%   framesPerEpochBin       = median(framesPerEpochBin(:), 'omitnan');
%   epochBinNoise           = [score.roi(activeCells).noise] / sqrt(framesPerEpochBin);
%   
%   %% Compute trial type (choice) selectivity 
%   [selectivity, preferredSide, ridgeToBkg]                        ...
%                           = computeChoiceSelectivity(epochDFF, sortTrials, cfg);
%     
%   %% Define reliability of firing as the fraction of trials with significant activity at the location of the peak response
%   [reliability, reliabilityCI, peakFrame, firingFieldRange]       ...
%                           = computeReliability(preferredSide, sortTrials, testTrials, epochDFF, epochBinNoise, cfg);
%   
%   %% Sort cells by location of peak activity within a given trial type
%   epochCenters            = toBinCenters([cfg.epochEdges, ceil(cfg.epochEdges(end))]);
%   peakEpoch               = epochCenters(peakFrame);
%   cellOrder               = arrayfun(@(x) find(preferredSide == x), enumeration('Choice'), 'UniformOutput', false)';
%   for iType = 1:numel(cellOrder)
%     iCell                 = cellOrder{iType};
%     [~,iOrder]            = sort(peakFrame(iCell));
%     cellOrder{iType}      = cellOrder{iType}(iOrder);
%   end
% 
% 
%   %%  ||    Shuffle test for significance of activity localization and stereotypy
%   %   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   
%   %% Decide how to parallelize tasks to minimize data transfer
%   pool                    = startParallelPool();
%   numChunks               = pool.NumWorkers;
%   numPerChunk             = round(linspace(1, cfg.numShuffles + 1, numChunks + 1));
%   numPerChunk(1)          = 1;
%   numPerChunk(end)        = cfg.numShuffles + 1;
%   numPerChunk             = diff(numPerChunk);
%   assert( sum(numPerChunk) == cfg.numShuffles );
%   
%   %% Setup random number generator and progress bar
%   randomGen               = RandStream.create('mlfg6331_64', 'Seed', 18217, 'NumStreams', numChunks);
%   RandStream.setGlobalStream(randomGen);
%   
%   %% Select only correct trials for the following, because inconsistent responses in error trials can affect metrics based on rotating data
% %   correctTrials           = sortTrials{end};
%   correctTrials           = [trialSet{:,1}];
%   isCorrect               = ismember(dataTrial, correctTrials);
%   correctDFF              = cellDFF(isCorrect,:);
%   correctBins             = dataBins(isCorrect,:);
%   
%   %% Define shuffle test using the same rotation factor for all cells (makes no difference because here cells are tested independently)
%   % Note that "reliability" doesn't make sense in the shuffled data because firing fields can be
%   % arbitrarily long, and therefore have high reliability just because it includes every bit of
%   % activity the cell might ever have. So only assess cells based on their ridge-to-background ratio
%   
%   shuffledRidge           = cell(numChunks, 1);
% %   shuffledReliability     = cell(numChunks, 1);
%   parfor iChunk = 1:numChunks
%     %% Generate random rotations using a reproducible random number stream
%     randStream            = RandStream.getGlobalStream();
%     randStream.Substream  = iChunk;
%     startFrame            = randi(size(correctDFF,1), 1, numPerChunk(iChunk));
%     chunkRidge            = nan(numPerChunk(iChunk),numel(ridgeToBkg));
% %     chunkConsist          = nan([numPerChunk(iChunk),size(reliability)]);
%     
%     %% Compute test statistics for shuffled experiments
%     for iExp = 1:numPerChunk(iChunk)
%       %% Generate null hypothesis pseudo-data by randomly rotating dF/F signals per cell
%       shuffledDFF         = correctDFF([startFrame(iExp):end,1:startFrame(iExp)-1], :);
%       shuffledDFF         = averageInBin(shuffledDFF, correctBins, nDataBins, 1);
% 
%       %% Compute ridge-to-background excess and reliability measures
%       [~, shuffledPref, chunkRidge(iExp,:)]        ...
%                                 = computeChoiceSelectivity(shuffledDFF, sortTrials, cfg);
% %       chunkConsist(iExp,:,:)    = computeReliability(shuffledPref, sortTrials, testTrials, shuffledDFF, epochBinNoise, cfg);
%     end
%     shuffledRidge{iChunk}       = chunkRidge;
% %     shuffledReliability{iChunk} = chunkConsist;
%   end
%   
%   shuffledRidge           = cat(1, shuffledRidge{:});
% %   shuffledReliability     = cat(1, shuffledReliability{:});
%   
%   %% Define significance as the fraction of shuffles with larger ridge/reliability metrics
%   ridgeSignificance       = mean(bsxfun(@gt, shuffledRidge, ridgeToBkg), 1);
% %   reliabilitySignificance = bsxfun(@gt, shuffledReliability, shiftdim(reliability,-1));
% %   reliabilitySignificance = squeeze(mean(reliabilitySignificance, 1));
% 
%   % Trial-vs-cell plotting code for testing purposes:
%   %{
%   [rr,ii]=sort(ridgeSignificance); tt=arrayfun(@(x,y,z,s,rs) sprintf('P(ridge)=%.3g, reliability=%.3g (P=%.3g), %s peak @ %d',x,y,rs,char(s),z), rr, reliability(ii), peakFrame(ii), preferredSide(ii), reliabilitySignificance(ii,end)', 'uniformoutput', false);
%   mm=MovieSlider(epochDFF(:,:,ii)); mm.setTitle(tt)
%   %}
% 
%   
%   %%  ||    Store output
%   %   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   
%   if exist(cfg.outputFile, 'file') == 2
%     movefile(cfg.outputFile, [cfg.outputFile '.old']);
%   end
%   
%   save( cfg.outputFile, 'versionInfo', 'globalID', 'trialIDs', 'trialSet', 'sortTrials', 'testTrials', 'epochDuration', 'epochDFF', 'basisKernel' ...
%       , 'selectivity', 'preferredSide', 'epochCenters', 'peakFrame', 'peakEpoch', 'firingFieldRange', 'cellOrder', 'ridgeToBkg', 'basisKernel2'   ...
%       , 'reliability', 'reliabilityCI', 'epochBinNoise', 'shuffledRidge', 'ridgeSignificance'                                                     ...
%       );
%   save( cfg.outputFile, '-append', '-struct', 'behavior');
%   
end


%---------------------------------------------------------------------------------------------------
function [selectivity, preferredSide, ridgeToBkg] = computeChoiceSelectivity(epochDFF, sortTrials, cfg)
  
  %% Compute trial type (choice) selectivity 
  selectivity             = harveyBinaryTypeSelectivity(epochDFF, sortTrials(1:2), cfg.minActiveFrac, cfg.minActiveFrames);
  preferredSide           = sign(selectivity);
  preferredSide( preferredSide == 0 ) = inf;
  preferredSide           = Choice( (preferredSide + 3) / 2 );
  
  %% Compute ridge-to-background excess using preferred-side trials per cell
  meanDFF                 = accumfun(2, @(x) mean(epochDFF(cfg.epochSel,x,:), 2, 'omitnan'), sortTrials);
  prefTrials              = sub2ind(butfirst(size(meanDFF)), min(double(preferredSide), cfg.numSides+1), 1:size(meanDFF,3));
  meanDFF                 = meanDFF(:,prefTrials);
  ridgeToBkg              = max(meanDFF,[],1) - halfSampleMode(meanDFF);
  
end

%---------------------------------------------------------------------------------------------------
function [reliability, reliabilityCI, peakFrame, firingFieldRange, epochBinNoise]   ...
                            = computeReliability(preferredSide, sortTrials, testTrials, epochDFF, epochBinNoise, cfg)

  %% Loop over cells to compute reliability of firing within a defined field
  peakFrame                 = nan(size(preferredSide));
  reliability               = nan(numel(preferredSide), numel(cfg.noiseFactor));
  reliabilityCI             = nan(numel(preferredSide), numel(cfg.noiseFactor), 2);
  firingFieldRange          = cell(numel(preferredSide),1);
  for iCell = 1:numel(preferredSide)
    %% Use preferred-side trials to define location of peak response for epoch-binned data 
    trials                  = sortTrials{min(preferredSide(iCell), cfg.numSides+1)};
    cellDFF                 = epochDFF(cfg.epochSel,trials,iCell);
    meanDFF                 = mean(cellDFF, 2, 'omitnan');
    [~,peakFrame(iCell)]    = max(meanDFF, [], 1, 'omitnan');
    
    %% Pad data circularly so that we can account for fields that span the end of a trial and the beginning of the next
    iPeak                   = peakFrame(iCell) + numel(meanDFF);      % shift this to the 2nd replicate
    meanDFF                 = repmat(meanDFF, 3, 1);
%     figure; plot(meanDFF);    hold on; plot(iPeak,meanDFF(iPeak),'*')
    
    %% Define the "firing field" of this cell as all consecutive bins within a given height of the peak
    isBaseline              = meanDFF < cfg.minFieldHeight * meanDFF(iPeak);
    firstFrame              = find(isBaseline(1:iPeak-1), 1, 'last') + 1;
    if isempty(firstFrame)
      firstFrame            = 1;
    end
    lastFrame               = iPeak + find(isBaseline(iPeak+1:end), 1, 'first') - 1;
    if isempty(lastFrame)
      lastFrame             = numel(meanDFF);
    end
    assert(~any(isBaseline(firstFrame:lastFrame)));
    
    %% Readjust range of frames for the firing field to account for the circular padding
    fieldRange              = (firstFrame:lastFrame) - size(cellDFF,1);
    invalid                 = fieldRange < 1;
    fieldRange(invalid)     = fieldRange(invalid) + size(cellDFF,1);
    invalid                 = fieldRange > size(cellDFF,1);
    fieldRange(invalid)     = fieldRange(invalid) - size(cellDFF,1);
    firingFieldRange{iCell} = fieldRange;
    
    %% Define reliability as fraction of test (cross-validation) trials where the sum activity within the firing field is significantly above noise
    trials                  = testTrials{min(preferredSide(iCell), cfg.numSides+1), 1};
    cellDFF                 = epochDFF(cfg.epochSel,trials,iCell);
    
    fieldNoise              = sqrt(numel(fieldRange)) * epochBinNoise(iCell);
    trialActivity           = sum(cellDFF(fieldRange,:), 1);
    sel                     = ~isnan(trialActivity);
    isConsistent            = trialActivity(sel)' > cfg.noiseFactor * fieldNoise;
    [reliability(iCell,:), reliabilityCI(iCell,:,:)]      ...
                            = binointerval(sum(isConsistent,1), sum(sel));
  end
  
end

