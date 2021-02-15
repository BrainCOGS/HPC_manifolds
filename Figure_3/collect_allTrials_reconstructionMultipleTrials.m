
%% Collects results

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

%% Setup
% Only animal E43
animal = 'E43';
fname = fnameStruct(3).fname;
load(fname);
deltaT = score.deltaT;

input_rng = 1;
rng(input_rng);

dimEmbed  = 5;

Ncrossval = 10;              % number of leave-out trials for crossvalidation
nleafs    = 500;             % list of leafs for fitting manifolds of different complexity
lmf       = 1;               % list of landmark_fractions

storage_directory  = 'D:\CrossValidation\';


%% Main Code

reconstructedAll = [];
sampleAll = [];
savepath = [storage_directory, animal '/'];

trialsList = [145 146 147 148 149 151 152 153 154 157 158 159 160];
folder_c = [savepath, 'reconstructContinuous/'];

for trial_idx = 1:length(trialsList)
    
    leaveout_trial = trialsList(trial_idx);
    folder_t = [folder_c, 'trialout_', num2str(leaveout_trial), '/'];
    
    f = [folder_t, 'minLeaves_' num2str(nleafs), '_lmf_', num2str(lmf), '_crossval.mat'];
    f2 = [folder_t, 'minLeaves_' num2str(nleafs), '_lmf_', num2str(lmf), '_manifold.mat'];
    disp(f);
    
    if isfile(f)
        load(f);
        
        % Set the number to 1 less than the dimensions desired
        reconstructedAll = [reconstructedAll; output_crossVal.reconstructedSave{dimEmbed-1}];
        lengthReconstructed(trial_idx) = size(output_crossVal.reconstructedSave{dimEmbed-1},1);
        
        sample = output_crossVal.sample;
        sampleAll = [sampleAll; sample];
    end
    
    if trial_idx==1
        load(f2);
    end
end


%% Collect and plot

% Get rid of neurons with no activity in window, doesn't matter much
% because only plotting 1-40, but nothing that low indices are removed
nic_output = extractVariables('all', 2, [145 146 147 148 149 151 152 153 154 157 158 159 160], 2, 0, 0, 5, [0 0], fname,'none','towers',1,1);
ROIactivities = nic_output.ROIactivities;
ROIactivities = ROIactivities(:,outMind.Neurons);

nic_output2 = extractVariables('all', 2, [145 146 147 148 149 151 152 153 154 157 158 159 160], 2, 0, 0, 5, [11 4], fname,'none','towers',1,1);
ROIactivities_11_4 = nic_output2.ROIactivities;
ROIactivities_11_4 = ROIactivities_11_4(:,outMind.Neurons);

% Get rid of timepoints that are all 0 across all neurons
ROIactivities_11_4_Datarange    = sum(ROIactivities_11_4,2)>0;
ROIactivities_11_4 = ROIactivities_11_4(ROIactivities_11_4_Datarange,:);
ROIactivities = ROIactivities(ROIactivities_11_4_Datarange,:);

% mind_preprocess the reconstructedAll and subtract the mean
reconstructedAll_11_4 = reconstructedAll - mean(reconstructedAll);
reconstructedAll_11_4 = mind_preprocess(reconstructedAll_11_4,11,4);


%% Plotting, the normalization is within the window

figure;
ax1 = subplot(2,2,1);
imagesc(mat2gray(ROIactivities(400:600,1:40)'));
xticks([0 5*(1/deltaT) 10*(1/deltaT)])
xticklabels({'0' '5' '10'});
set(gca, 'box', 'off')
ylabel('Neuron #');
xlabel('Time (s)');
title('raw')
axis square

ax2 = subplot(2,2,2);
imagesc(mat2gray(ROIactivities_11_4(400:600,1:40)'));
xticks([0 5*(1/deltaT) 10*(1/deltaT)])
xticklabels({'0' '5' '10'});
set(gca, 'box', 'off')
ylabel('Neuron #');
xlabel('Time (s)');
title('raw (smoothed and thresholded)')
axis square

ax3 = subplot(2,2,3);
imagesc(mat2gray(reconstructedAll(400:600,1:40)'));
xticks([0 5*(1/deltaT) 10*(1/deltaT)])
xticklabels({'0' '5' '10'});
set(gca, 'box', 'off')
xlabel('Time (s)');
title('reconstructed')
axis square

ax4 = subplot(2,2,4);
imagesc(mat2gray(reconstructedAll_11_4(400:600,1:40)'));
xticks([0 5*(1/deltaT) 10*(1/deltaT)])
xticklabels({'0' '5' '10'});
set(gca, 'box', 'off')
xlabel('Time (s)');
title('reconstructed (smoothed and thresholded)')
axis square

Link = linkprop([ax1 ax2 ax3 ax4],{'XLim', 'YLim', 'ZLim'});

