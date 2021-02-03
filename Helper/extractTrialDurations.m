function outputTrialDurations = extractTrialDurations(fn_modeling)
% fn_modeling = 'C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat'
load(fn_modeling)
keepTrials = find([score.trial.mainTrial]==1 & [score.trial.goodQuality]==1);
trial_keepTrials = trial(keepTrials);

for i=1:length(keepTrials)
   
    tTrialStart2tTrialEnd(i) = trial_keepTrials(i).tTrialEnd - trial_keepTrials(i).tTrialStart;
    tCueEntry2tMemEntry(i)   = trial_keepTrials(i).tMemEntry - trial_keepTrials(i).tCueEntry;
    tMemEntry2tArmEntry(i)   = trial_keepTrials(i).tArmEntry - trial_keepTrials(i).tMemEntry;
    tTrialStart2tCueEntry(i) = trial_keepTrials(i).tCueEntry - trial_keepTrials(i).tTrialStart;
    tArmEntry2tTrialEnd(i)   = trial_keepTrials(i).tTrialEnd - trial_keepTrials(i).tArmEntry; 
    
end

outputTrialDurations.tTrialStart2tTrialEnd = tTrialStart2tTrialEnd;
outputTrialDurations.tCueEntry2tMemEntry = tCueEntry2tMemEntry;
outputTrialDurations.tMemEntry2tArmEntry = tMemEntry2tArmEntry;
outputTrialDurations.tTrialStart2tCueEntry = tTrialStart2tCueEntry;
outputTrialDurations.tArmEntry2tTrialEnd = tArmEntry2tTrialEnd;

outputTrialDurations.meanLength = mean(tTrialStart2tTrialEnd);
outputTrialDurations.meanCue    = mean(tCueEntry2tMemEntry);
outputTrialDurations.meanMem    = mean(tMemEntry2tArmEntry);
outputTrialDurations.meanStart  = mean(tTrialStart2tCueEntry);
outputTrialDurations.meanArm    = mean(tArmEntry2tTrialEnd);

