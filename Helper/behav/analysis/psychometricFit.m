function psych = psychometricFit(choice,nCues_RminusL,fitNotBinned,bins,fitBinned)

% psych = psychometricFit(choice,nCues_RminusL,fitNotBinned,bins,fitBinned)

%% defaults
if nargin < 3
  fitNotBinned = 1;
end
if nargin < 4 || isempty(bins)
  bins = analysisParams.psychbins;
end
if nargin < 5
  fitBinned = 1;
end

%% psychometric curve from raw data (all diffs and binned)
psych.perfPsych_bins      = bins;
psych.stdInterval         = normcdf(1, 0, 1) - normcdf(-1, 0, 1);
psych.perfPsych_xaxis     = -15:15;

psych.perfPsych           = zeros(size(psych.perfPsych_xaxis));
psych.nR                  = psych.perfPsych;
psych.nTrials             = psych.perfPsych;
psych.perfPsych_xaxisBins = zeros(1,numel(psych.perfPsych_bins)-1);
psych.perfPsych_binned    = zeros(size(psych.perfPsych_xaxisBins));
psych.nR_binned           = psych.perfPsych_binned;
psych.nTrials_binned      = psych.perfPsych_binned;

% all #r - #l, fraction went right
for ii = 1:length(psych.perfPsych_xaxis)
  psych.nR(ii)        = sum(choice(nCues_RminusL==psych.perfPsych_xaxis(ii)) == analysisParams.rightCode);
  psych.nTrials(ii)   = sum(nCues_RminusL==psych.perfPsych_xaxis(ii));
  psych.perfPsych(ii) = psych.nR(ii)/psych.nTrials(ii);
end
% binned
for ii = 1:length(psych.perfPsych_bins)-1
  psych.nR_binned(ii)           = sum(choice(nCues_RminusL>=psych.perfPsych_bins(ii) & ...
    nCues_RminusL<psych.perfPsych_bins(ii+1)) == analysisParams.rightCode);
  psych.nTrials_binned(ii)      = sum(nCues_RminusL>=psych.perfPsych_bins(ii)        & ...
    nCues_RminusL<psych.perfPsych_bins(ii+1));
  psych.perfPsych_binned(ii)    =  psych.nR_binned(ii)/psych.nTrials_binned(ii);
  
  % bin center needs to be weighted average of delta per ntrials
  dvals = psych.perfPsych_bins(ii):psych.perfPsych_bins(ii+1) -1;
  for jj = 1:numel(dvals); nt(jj) = sum(nCues_RminusL==dvals(jj)); end
  try
    psych.perfPsych_xaxisBins(ii) = sum(dvals.*nt)./sum(nt);
  catch
    psych.perfPsych_xaxisBins(ii) = psych.perfPsych_bins(ii) + mode(diff(psych.perfPsych_bins))/2;
  end
end

% psychometric curve using jeffrey's method (all diffs and binned)
% [psych.perfPsychJ,psych.perfPsychJSEM]               = ...
%     binofit(psych.nR,psych.nTrials,1-psych.stdInterval);
% [psych.perfPsychJ_binned,psych.perfPsychJ_binnedSEM] = ...
%     binofit(psych.nR_binned,psych.nTrials_binned,1-psych.stdInterval);
[psych.perfPsychJ,psych.perfPsychJSEM]               = ...
  binointervalLP(psych.nR,psych.nTrials,1-psych.stdInterval);
[psych.perfPsychJ_binned,psych.perfPsychJ_binnedSEM] = ...
  binointervalLP(psych.nR_binned,psych.nTrials_binned,1-psych.stdInterval);

% 4-parameter sigmoid curve
% Sigmoid function fit, weighted by 1/sigma^2 where sigma is the symmetrized error
psych.sigmoid         = @(O,A,lambda,x0,x) O + A ./ (1 + exp(-(x-x0)/lambda));
psych.sigmoidSlope    = @(A,lambda) A ./ (4*lambda);

% fit on binned
if fitBinned
  try
    psych.fit.xaxis       = psych.perfPsych_bins(1):0.05:psych.perfPsych_bins(end);
    psych.fit.params      = fit ( psych.perfPsych_xaxisBins(~isnan(psych.perfPsych_binned))', ...
      psych.perfPsych_binned(~isnan(psych.perfPsych_binned))',    ...
      psych.sigmoid                                               ...
      , 'StartPoint'      , [0 1 8 0]                             ...
      , 'Weights'         , ((psych.perfPsychJ_binnedSEM(~isnan(psych.perfPsych_binned),2)- psych.perfPsychJ_binnedSEM(~isnan(psych.perfPsych_binned),1)) / 2).^-2   ...
      , 'MaxIter'         , 400 ...
      );
    %                         , 'Lower'           , [-0.5 -2 0  data.deltaRLBins(1)  ]       ...
    %                         , 'Upper'           , [ 0.5  2 10 data.deltaRLBins(end)]       ...
    
    psych.fit.stdInt      = predint(psych.fit.params, 0, psych.stdInterval, 'functional');
    psych.fit.bias        = psych.fit.params(0);
    psych.fit.biasCI      = predint(psych.fit.params, 0, 0.95, 'functional')';
    psych.fit.biasErr     = (psych.fit.stdInt(2) - psych.fit.stdInt(1)) / 2;
    psych.fit.slope       = psych.sigmoidSlope(psych.fit.params.A, psych.fit.params.lambda);
    psych.fit.curve       = feval(psych.fit.params,psych.fit.xaxis');
  catch
    psych.fit = [];
  end
end

% fit on full
if fitNotBinned
  try
    psych.fitAll.xaxis       = psych.perfPsych_xaxis(1):0.05:psych.perfPsych_xaxis(end);
    psych.fitAll.params      = fit ( psych.perfPsych_xaxis(~isnan(psych.perfPsychJ))', ...
      psych.perfPsychJ(~isnan(psych.perfPsychJ)),    ...
      psych.sigmoid                                               ...
      , 'StartPoint'      , [0 1 8 0]                             ...
      , 'Weights'         , ((psych.perfPsychJSEM(~isnan(psych.perfPsychJ),2)- psych.perfPsychJSEM(~isnan(psych.perfPsychJ),1)) / 2).^-2   ...
      );
    %                         , 'Lower'           , [-0.5 -2 0  data.deltaRLBins(1)  ]       ...
    %                         , 'Upper'           , [ 0.5  2 10 data.deltaRLBins(end)]       ...
    
    psych.fitAll.stdInt      = predint(psych.fitAll.params, 0, psych.stdInterval, 'functional');
    psych.fitAll.bias        = psych.fitAll.params(0);
    psych.fitAll.biasCI      = predint(psych.fitAll.params, 0, 0.95, 'functional')';
    psych.fitAll.biasErr     = (psych.fitAll.stdInt(2) - psych.fitAll.stdInt(1)) / 2;
    psych.fitAll.slope       = psych.sigmoidSlope(psych.fitAll.params.A, psych.fitAll.params.lambda);
    psych.fitAll.curve       = feval(psych.fitAll.params,psych.fitAll.xaxis');
  catch
    psych.fitAll = [];
  end
end
