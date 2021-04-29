
%% Generate the example that's identical method to the reconstruction in Figure 3

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

outputReconstructionMultipleROIs_E43_1to40 = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, '\\192.168.0.233\Neuroscience\CrossValidation_holdOneCell\', 'E43', [1:40]);

sourceData_S4d = plot_allTrials_reconstructionMultipleROIs_Function(outputReconstructionMultipleROIs_E43_1to40, [145 146 147 148 149 151 152 153 154 157 158 159 160], fnameStruct(3).fname);


%% Get the correlation coefficient in the 10 held-out ROIs in all animals

% Used the .sh files, i.e. "M:\enieh\mind\mind_collect_reconstructionROIs_E43_spock.sh" to
% generate the outputReconstructionMultipleROIs files
% Needed to manually add the 10 trials for E43, since they're chosen
% randomly from the top 25 but also in first 40 (this is taken care of in
% the .sh file already

% OLD ****************************************************
% Editted code to fix the idim starting at 1, but this dataset has the 
% additional column, so I manually editted out the extra column and saved
% as the data from the folder: C:\Neuroscience\imaging\FINAL\reconstructROIs_Data
% outputReconstructionMultipleROIs_E43 = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, '\\192.168.0.233\Neuroscience\CrossValidation_holdOneCell\', 'E43', [2 6 8 9 14 20 26 28 31 34]);
% for i=1:length(fnameStruct)
%     if i==3
%         outputReconstructionMultipleROIs{i} = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, '\\192.168.0.233\Neuroscience\CrossValidation_holdOneCell\', i, [2 6 8 9 14 20 26 28 31 34]);
%     else
%         outputReconstructionMultipleROIs{i} = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, '\\192.168.0.233\Neuroscience\CrossValidation_holdOneCell\', i, []);
%     end
%     disp(['Finished animal ' num2str(i) 'of ' num2str(length(fnameStruct))]);
%     save('reconstructMultipleROIs_intermediateSave_TRASH.mat', '-v7.3');
%     % NEED TO FIX TO MANUALLY SET ROIs FOR E43
%     %***************************************************
% end
% OLD END *************************************************

load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E22_spock.mat')
load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E39_spock.mat')
load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E43_spock.mat')
load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E44_spock.mat')
load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E47_spock.mat')
load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E48_spock.mat')
load('C:\Neuroscience\imaging\FINAL\reconstructROIs_Data\outputReconstructionMultipleROIs_E65_spock.mat')

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

corrAll(1,:) = mean(outputReconstructionMultipleROIs_E22.corrAll,1);
corrAll(2,:) = mean(outputReconstructionMultipleROIs_E39.corrAll,1);
corrAll(3,:) = mean(outputReconstructionMultipleROIs_E43.corrAll,1);
corrAll(4,:) = mean(outputReconstructionMultipleROIs_E44.corrAll,1);
corrAll(5,:) = mean(outputReconstructionMultipleROIs_E47.corrAll,1);
corrAll(6,:) = mean(outputReconstructionMultipleROIs_E48.corrAll,1);
corrAll(7,:) = mean(outputReconstructionMultipleROIs_E65.corrAll,1);

figure;
nieh_barSEM(corrAll);
sourceData_S4e = corrAll';
hold on;
        scatter([ones(length(fnameStruct),1); ...
                 ones(length(fnameStruct),1)*2; ...
                 ones(length(fnameStruct),1)*3; ...
                 ones(length(fnameStruct),1)*4; ...
                 ones(length(fnameStruct),1)*5; ...
                 ones(length(fnameStruct),1)*6; ...
                 ],corrAll(:), '.');
ylabel('Mean Cross-Validated Score');
xticklabels({'2', '3', '4', '5', '6', '7'});
xlabel('Num Embedding Dims');
set(gca, 'box', 'off')
axis square
