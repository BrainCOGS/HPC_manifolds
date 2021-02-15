function outputMakeVAPlot = mind_makeVAplot(fnameStruct, taskType)

% Make the plot of View Angle as function of Position

%% Save inputs

outputMakeVAPlot.fnameStruct = fnameStruct;
outputMakeVAPlot.taskType    = taskType;

%% Do function for all animals

left_VA_all = [];
right_VA_all = [];
for i=1:length(fnameStruct)
   
    fname = fnameStruct(i).fname;
    
    [meanL(i,:), meanR(i,:), semL(i,:), semR(i,:), left_VA_cur, right_VA_cur] = calcVAmean(fname, taskType);
    left_VA_all = [left_VA_all left_VA_cur];
    right_VA_all = [right_VA_all right_VA_cur];

end

outputMakeVAPlot.meanL = meanL;
outputMakeVAPlot.meanR = meanR;
outputMakeVAPlot.semL  = semL;
outputMakeVAPlot.semR  = semR;
outputMakeVAPlot.left_VA_all = left_VA_all;
outputMakeVAPlot.right_VA_all = right_VA_all;

end


%% Load Variables depending on taskType

function [meanL, meanR, semL, semR, left_VA, right_VA] = calcVAmean(fname, taskType)

if strcmp(taskType, 'alternation')==1 || strcmp(taskType, 'Alternation')==1
    nic_output = extractVariables('all', 2, 'goodTrials', 2, 30, 0, 11, [11 4], fname,'none','alternation', 1, 1);
    
elseif strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 30, 0, 11, [11 4], fname,'none', 'towers', 1, 1);
    
elseif strcmp(taskType, 'alternationJeff')==1 || strcmp(taskType, 'AlternationJeff')==1
    nic_output = extractVariables('all', 6, 'goodTrials', 2, 30, 0, 11, [11 4], fname,'none','alternationJeff', 1, 1);
    
end

behavioralVariables = nic_output.behavioralVariables;

trial=behavioralVariables.Trial;
VA   =rad2deg(behavioralVariables.ViewAngle);
C    =behavioralVariables.Choice;

VA_shape = reshape(VA, [31 length(unique(trial))]);
C_shape = reshape(C, [31 length(unique(trial))]);

C_trial = C_shape(1,:);

right_VA = VA_shape(:,C_trial==1);
left_VA = VA_shape(:,C_trial==0);

meanL = mean(left_VA,2);
meanR = mean(right_VA,2);
semL  = nieh_sem(left_VA');
semR  = nieh_sem(right_VA');

end