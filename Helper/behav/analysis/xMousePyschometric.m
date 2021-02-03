function psych = xMousePyschometric(choice,nCues_RminusL,mouseID)

% psych = xMousePyschometric(choice,nCues_RminusL,mouseID)

psych.mice            = unique(mouseID);
psych.nmice           = length(psych.mice);
psych.slopes          = nan(1,psych.nmice);
psych.perfOv          = nan(1,psych.nmice);
psych.lapse           = nan(1,psych.nmice);
tt                    = zeros(size(choice));
tt(nCues_RminusL > 0) = 1;

tempPerf = []; tempfit = []; tempSEM = [];
for ii = 1:psych.nmice
    psych.mouse(ii)   = psychometricFit(choice(mouseID == psych.mice(ii)),nCues_RminusL(mouseID == psych.mice(ii)),1);
    psych.nTrials(ii) = sum(psych.mouse(ii).nTrials);
    if psych.nTrials(ii) > 500
        tempPerf(:,end+1)  = psych.mouse(ii).perfPsychJ;
        tempSEM(:,:,end+1) = psych.mouse(ii).perfPsychJSEM;
        tempfit(:,end+1)   = psych.mouse(ii).fitAll.curve;
        psych.slopes(ii)   = psych.mouse(ii).fitAll.slope;
        psych.perfOv(ii)   = 100.*(sum(choice(mouseID == psych.mice(ii))==tt(mouseID == psych.mice(ii)))/sum(mouseID == psych.mice(ii)));
        psych.lapse(ii)    = 100 - 100.*(sum(choice(mouseID == psych.mice(ii) & abs(nCues_RminusL) >= 10)==tt(mouseID == psych.mice(ii) & abs(nCues_RminusL) >= 10))./sum(mouseID == psych.mice(ii) & abs(nCues_RminusL) >= 10));
    end
end
psych.goodCurves = tempPerf;
psych.goodFits   = tempfit;
psych.goodSEM    = tempSEM;
psych.fitXaxis   = psych.mouse(ii).fitAll.xaxis;
psych.xmouseMean = nanmean(tempPerf,2);
psych.xmouseSEM  = nanstd(tempPerf,0,2)./sqrt(psych.nmice-1);

