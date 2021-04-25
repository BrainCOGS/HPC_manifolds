function sourceData_S4d = plot_allTrials_reconstructionMultipleROIs_Function(outputReconstructionMultipleROIs, trialsList, fname)
% trialsList = [145 146 147 148 149 151 152 153 154 157 158 159 160];

load(fname);
deltaT = score.deltaT;

reconstructedAll = outputReconstructionMultipleROIs.reconstructedAll;
corrAll = outputReconstructionMultipleROIs.corrAll;

nic_outputFull = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [0 0], fname,'none','towers', 1, 1);
behavioralVariablesFull = nic_outputFull.behavioralVariables;
ROIactivitiesFull = nic_outputFull.ROIactivities;

nic_output = extractVariables('all', 2, trialsList, 2, 0, 0, 5, [0 0], fname,'none','towers',1,1);
ROIactivities = nic_output.ROIactivities;

nic_output2 = extractVariables('all', 2, trialsList, 2, 0, 0, 5, [11 4], fname,'none','towers',1,1);
ROIactivities_11_4 = nic_output2.ROIactivities;

ROIactivities_11_4_Datarange    = sum(ROIactivities_11_4,2)>0;
ROIactivities_11_4 = ROIactivities_11_4(ROIactivities_11_4_Datarange,:);
ROIactivities = ROIactivities(ROIactivities_11_4_Datarange,:);

% Get the reconstructedAll data that's only in the desired trials
trialsKeep = ismember(behavioralVariablesFull.Trial,trialsList);
justTrialData = reconstructedAll(trialsKeep,:);
justTrialData = justTrialData(ROIactivities_11_4_Datarange,:);
reconstructedAll_11_4 = reconstructedAll - mean(reconstructedAll);
reconstructedAll_11_4 = mind_preprocess(reconstructedAll_11_4,11,4);
justTrialData_11_4 = reconstructedAll_11_4(trialsKeep,:);
justTrialData_11_4 = justTrialData_11_4(ROIactivities_11_4_Datarange,:);

mean(corrAll)
std(corrAll)

figure;
ax1 = subplot(1,2,1);
imagesc(mat2gray(ROIactivities_11_4(400:600,1:40)'));
xticks([0 5*(1/deltaT) 10*(1/deltaT)])
xticklabels({'0' '5' '10'});
set(gca, 'box', 'off')
ylabel('Neuron #');
xlabel('Time (s)');
title('raw (smoothed and thresholded)')
axis square

ax2 = subplot(1,2,2);
imagesc(mat2gray(justTrialData_11_4(400:600,1:40)'));
xticks([0 5*(1/deltaT) 10*(1/deltaT)])
xticklabels({'0' '5' '10'});
set(gca, 'box', 'off')
xlabel('Time (s)');
title('reconstructed (smoothed and thresholded)')
axis square

sourceData_S4d = mat2gray(justTrialData_11_4(400:600,1:40)');