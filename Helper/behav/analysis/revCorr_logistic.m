function rc = revCorr_logistic(trialtype,choice,cuePos_R,cuePos_L,info,subtractExpected,nbins,doBoot,nBootIter)

% rc = revCorr_logistic(trialtype,choice,cuePos_R,cuePos_L,info,subtractExpected,nbins)
% calculates reverse correlation using logistic regression, #R-#L/spatial
% bin as predictors
%
% LP feb 2016

warning('off','all'); 

%% defaults
if nargin < 5; info             = [];      end
if nargin < 6; subtractExpected = 0;       end
if nargin < 7; nbins            = 5;       end
if nargin < 8; doBoot           = false;   end
if nargin < 9; nBootIter        = 200;     end

%% cue prob, avg number of cues etc
if isempty(info);
  rc.lCue = 200;
else
  try
    rc = getCueStats(info,cuePos_R,cuePos_L);
  catch
    rc = getCueStatsFromMaze(info,cuePos_R,cuePos_L); % in case info is actually already the maze structure
  end
end

%% make R - L predictor matrix, trials x spatial bins
rc.bins = linspace(10,rc.lCue,nbins+1);
ntrials = numel(choice);
RminusLmat = zeros(ntrials,nbins);
for ii = 1:ntrials
  rhist = histc(cuePos_R{ii},rc.bins);
  lhist = histc(cuePos_L{ii},rc.bins);
  RminusLmat(ii,:) = rhist(1:end-1) - lhist(1:end-1);
end

% delete unfinished trials
delidx = find(choice~=analysisParams.rightCode & choice~=analysisParams.leftCode);
RminusLmat(delidx,:) = [];
choice(delidx)       = [];
trialtype(delidx)    = [];

% subtracted expected number of cues
if subtractExpected
  rc.meanNumCues(2)  = rc.cueMeanCt/(1+exp(rc.cueProb));
  rc.meanNumCues(1)  = rc.cueMeanCt-rc.meanNumCues(2);
  rc.expectedNumCues = -sign(trialtype-.5)*(rc.meanNumCues(1)-rc.meanNumCues(2))/(length(rc.bins)-1);
  predMat            = RminusLmat-repmat(rc.expectedNumCues',1,length(rc.bins)-1);
else
  predMat = RminusLmat;
end

%% fit logistic regression using glm fit
% choice(choice==0)        = -1;
[rc.weights,~, rc.stats] = glmfit(double(predMat),double(choice'),'binomial','link','logit');
rc.values                = rc.weights(2:end);
rc.decayIndex            = mean(rc.values(end-1:end))/mean(rc.values(1:2));

%% do bootstrap for error bars if necessary
if doBoot
  values  = zeros(nBootIter,numel(rc.values));
  ind     = zeros(nBootIter,1);
  ntrials = numel(choice);
  rng('default')
  for iBoot = 1:nBootIter
    idx     = randsample(ntrials,ntrials,true); % sample with replacement from all trials
    weights = glmfit(double(predMat(idx,:)),double(choice(idx)'),'binomial','link','logit');
    values(iBoot,:) = weights(2:end);
    
    mat = predMat;
    for iTrial = 1:ntrials; mat(iTrial,:) = predMat(iTrial,randperm(nbins)); end
    weights = glmfit(double(mat),double(choice'),'binomial','link','logit');
    weights = weights(2:end);
    ind(iBoot) = mean(weights(end-1:end))/mean(weights(1:2));
  end
  rc.sem          = std(values);
  if rc.decayIndex <= mean(ind)
    rc.decayIndex_p = sum(ind<rc.decayIndex)/nBootIter;
  else
    rc.decayIndex_p = sum(ind>rc.decayIndex)/nBootIter;
  end
  
else
  rc.sem          = rc.stats.se(2:end);
  rc.decayIndex_p = nan;
end

warning('on','all'); 

end

%%
function rc = getCueStats(info,cuePos_R,cuePos_L)

mazes        = feval(info(1).protocol);
rc.lCue      = mazes(info.mainMazeID).lCue;
rc.cueProb   = mazes(info.mainMazeID).cueProbability;

ntrials      = length(cuePos_R);
numcues      = zeros(1,ntrials);
for ii = 1:ntrials
  numcues(ii) = numel(cuePos_R{ii}) + numel(cuePos_L{ii});
end
rc.cueMeanCt = mean(numcues);

end

%%
function rc = getCueStatsFromMaze(mazes,cuePos_R,cuePos_L)

rc.lCue      = mazes.mazes(mazes.mainMazeID).lCue;
rc.cueProb   = mazes.mazes(mazes.mainMazeID).cueProbability;

ntrials      = length(cuePos_R);
numcues      = zeros(1,ntrials);
for ii = 1:ntrials
  numcues(ii) = numel(cuePos_R{ii}) + numel(cuePos_L{ii});
end
rc.cueMeanCt = mean(numcues);

end