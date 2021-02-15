function [triallength, err] = TrialLengthFromD(fname, taskType, outMind)
% taskType can be 'towers', 'alternation', or 'alternationJeff'

D = outMind.dat.forestdat.rwd.Dg;
landmarks = outMind.dat.forestdat.lm.idx;

if strcmp(taskType,'towers')==1 || strcmp(taskType, 'T7')==1 || strcmp(taskType,'Towers')==1
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);
elseif strcmp(taskType, 'alternation')==1 || strcmp(taskType, 'Alternation')==1
    nic_output = extractVariables('all', 2, 'goodTrials', 2, 0, 0, 11, [11 4], fname,'none','alternation', 1, 1);
elseif strcmp(taskType, 'alternationJeff')==1
    nic_output = extractVariables('all', 6, 'goodTrials', 2, 0, 0, 11, [11 4], fname,'none','alternation', 1, 1);
end

behavioralVariables = nic_output.behavioralVariables;
behavioralVariables = behavioralVariables(outMind.Datarange,:);

trialn = behavioralVariables.Trial;

trial_of_landmarks         = trialn(landmarks);
utrials = unique(trial_of_landmarks);

ll = zeros(length(unique(trial_of_landmarks)), 1);
for t_idx1 = 1:length(utrials)
    tt1 = utrials(t_idx1);
    entries1 = trial_of_landmarks == tt1;
    ll(t_idx1) = D( min(find(entries1)), max(find(entries1)));
end

triallength = mean(ll);
err = std(ll)/sqrt(length(ll));
