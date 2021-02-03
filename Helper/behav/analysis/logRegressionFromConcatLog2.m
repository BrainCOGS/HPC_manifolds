function rc = logRegressionFromConcatLog2(lg,minNT,nBins)

% rc = logRegressionFromConcatLog(lg,minNT)

if nargin < 2
  minNT    = 1000;
end
if nargin < 3
  nBins    = 5;
end

% lg         = cleanupConcatLog(lg,minNT);
mice       = unique(lg.mouseID);
nmice      = numel(mice);

rc.overall = revCorr_logistic(lg.trialType,lg.choice,lg.cuePos_R,lg.cuePos_L,[],0,nBins);
rc.bins    = rc.overall.bins;
rc.values  = rc.overall.values;

% stats per mouse
% rc.goodMiceID            = analysisParams.miceBehav(mice);
rc.ntrials               = zeros(1,nmice);
rc.goodValues            = zeros(numel(rc.values),nmice);
rc.slopes                = zeros(1,nmice);
rc.intercepts            = zeros(1,nmice);
rc.slopesPoly1           = zeros(1,nmice);
rc.interceptsPoly1       = zeros(1,nmice);
rc.quadratic             = zeros(1,nmice);
rc.peaks                 = zeros(1,nmice);
rc.isSigSlope            = zeros(1,nmice);
rc.isSigSlopeBoot        = zeros(1,nmice);
rc.isSigQuadr            = zeros(1,nmice);
rc.nSessGoodMice         = zeros(1,nmice);
rc.r2lin                 = zeros(1,nmice);
rc.r2exp                 = zeros(1,nmice);
rc.expFactor             = zeros(1,nmice);
rc.derivQuadr            = zeros(1,nmice);
rc.decayIndex            = zeros(1,nmice);
rc.isSigDerivQuadrBoot   = zeros(1,nmice);
rc.isSigDecayIndexBoot   = zeros(1,nmice);


% for individual mice
for mm = 1:numel(mice)
  [choice, cpr, cpl, tt]    ...
                 = selectMouseTrials(lg, mice(mm), 'choice', 'cuePos_R', 'cuePos_L', 'trialType');
  rc.mice(mm)    = revCorr_logistic(tt, choice, cpr, cpl ,[],0,nBins,true); 
  rc.ntrials(mm) = numel(choice);
  
  % general info
  rc.nSessGoodMice(:,mm) = numel(unique(lg.sessionID(lg.mouseID == mice(mm))));
  x                      = rc.bins(1:end-1)'+mode(diff(rc.bins))/2;
  y                      = rc.mice(mm).values;
  rc.decayIndex(:,mm)    = rc.mice(mm).decayIndex;%nanmean(y(end-1:end))/nanmean(y(1:2));
  rc.isSigDecayIndexBoot(:,mm) = rc.mice(mm).decayIndex_p<.05;
  rc.goodValues(:,mm)    = rc.mice(mm).values;
  rc.peaks(:,mm)         = x(y==max(y));
  
  % quadratic fit
  [temp,stats]            = fit(x,y,'poly2');
  yhat                    = fiteval(temp,x);
  rc.derivQuadr(:,mm)     = nanmean(diff(yhat));
  coeffs                  = coeffvalues(temp);
  ci                      = confint(temp);
  rc.r2poly2(mm)          = stats.rsquare;
  rc.fitStats_poly2(mm)   = stats;
  rc.slopes(:,mm)         = coeffs(2);
  rc.intercepts(:,mm)     = coeffs(3);
  rc.quadratic(:,mm)      = coeffs(1);
  if ci(1,2) <= 0 && ci(2,2) >= 0; sig = 0; else sig = 1; end
  if ci(1,1) <= 0 && ci(2,1) >= 0; sig2 = 0; else sig2 = 1; end
  rc.isSigSlope(:,mm)     = sig;
  rc.isSigQuadr(:,mm)     = sig2;
   
%   % bootstrap for slope significance
%   niter = 500;
%   sli   = zeros(1,niter);
%   deri  = zeros(1,niter);
%   indi  = zeros(1,niter);
%   sem   = zeros(niter,numel(rc.values));
%   for iBoot = 1:niter %shuffle bin ids
%     idx          = randperm(numel(y));%randsample(numel(tt),numel(tt),true);
%     irc          = y(idx); %revCorr_logistic(tt(idx),choice(idx),cpr(idx),cpl(idx),[],0,nBins);
%     indi(iBoot)  = nanmean(irc(end-1:end))./nanmean(irc(1:2));
%     temp         = fit(x,y,'poly2');%irc.bins(1:end-1)'+mode(diff(irc.bins))/2,irc.values,'poly2');
%     yh           = fiteval(temp, x);%irc.bins(1:end-1)'+mode(diff(irc.bins))/2);
%     deri(iBoot)  = nanmean(diff(yh));
%     cf           = coeffvalues(temp);
%     sli(iBoot)   = cf(2);
%     sem(iBoot,:) = irc;%.values;
%   end
%   tempStd(mm).bootStd = std(sem);
%   if sum(sign(sli)~=sign(coeffs(2))) > .05*niter; sig = 0; else sig = 1; end
%   rc.isSigSlopeBoot(:,mm)      = sig;
%   
%   if sum(sign(deri)~=sign(rc.derivQuadr(mm))) > .05*niter; sig = 0; else sig = 1; end
%   rc.isSigDerivQuadrBoot(:,mm) = sig;
%   
%   if sum(sign(indi)~=sign(rc.decayIndex(mm))) > .05*niter; sig = 0; else sig = 1; end
%   rc.isSigDecayIndexBoot(:,mm) = sig;
  
  % linear fit
  [temp,stats]                = fit(x,y,'poly1');
  coeffs                      = coeffvalues(temp);
  rc.r2lin(mm)                = stats.rsquare;
  rc.fitStats_lin(mm)         = stats;
  rc.slopesPoly1(:,mm)        = coeffs(1);
  rc.interceptsPoly1(:,mm)    = coeffs(2);
  
  % exponential fit
  [temp,stats]                = fit(rc.bins(1:end-1)'+mode(diff(rc.bins))/2,rc.mice(mm).values,...
    fittype(@(a,b,c,x) a+b*exp(c*x)),'startPoint',[.1 0 .2],'lower',[0 -5 -1],'upper',[1 5 1]);
  temp                        = coeffvalues(temp);
  rc.r2exp(mm)                = stats.rsquare;
  rc.fitStats_exp(mm)         = stats;
  rc.expFactor(:,mm)          = temp(3);
  %             xaxis = rc.bins(1:end-1)'+mode(diff(rc.bins))/2;
  %             figure; plot(xaxis,rc.mice(mm).values); hold on; plot(xaxis,rc.slopes(end).*xaxis+rc.intercepts(end)); pause; close;
end

rc.nGoodMice  = size(rc.goodValues,2);
%rc.sem        = nanstd(rc.goodValues,0,2)./sqrt(rc.nGoodMice-1);
rc.sem        = nieh_sem(rc.goodValues')';

for mm = 1:numel(rc.mice)
  rc.mice(mm).bootStd = rc.mice(mm).sem;
end


% plotRevCorr_ctrl(rc)

%inspect results
% mcount = 1;
% for ii = 1:numel(rc.mice)
%    if rc.ntrials(ii) > 1000
%       figure; hold on
%       errorbar(20:40:180,rc.mice(ii).values,rc.mice(ii).bootStd,'k-')
%       x = 0:10:200;
%       y = rc.intercepts(mcount) + rc.slopes(mcount).*x + rc.quadratic(mcount).*(x.^2);
%       plot(x,y,'r-'); ylim([0 .35])
%       title(sprintf('%s, deriv = %1.2g (sig = %d)\nindex = %1.2g (sig = %d)',rc.goodMiceID{mcount},rc.derivQuadr(mcount),rc.isSigDerivQuadrBoot(mcount),rc.decayIndex(mcount),rc.isSigDecayIndexBoot(mcount)))
%       mcount = mcount+1;
%       pause; close
%    end
% end

