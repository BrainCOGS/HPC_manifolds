
%% Generate the example that's identical to the reconstruction in Figure 3

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

outputReconstructionMultipleROIs_E43_1to40 = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, 'D:\CrossValidation_holdOneCell\', 'E43', [1:40]);


%% Get the correlation coefficient in the 10 held-out ROIs in all animals

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

for i=1:length(fnameStruct)
    outputReconstructionMultipleROIs{i} = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, 'D:\CrossValidation_holdOneCell\', i, []);
    disp(['Finished animal ' num2str(i) 'of ' num2str(length(fnameStruct))]);
    % NEED TO FIX TO MANUALLY SET ROIs FOR E43
    %***************************************************
end

% Used the .sh files, i.e. mind_collect_reconstructionROIs_E43.sh
% Needed to manually add the 10 trials for E43, since they're chosen
% randomly from the top 25 but also in first 40

% Load the data from the folder: C:\Neuroscience\imaging\FINAL\reconstructROIs_Data

%% Editted code to fix the idim starting at 1, manually editted data
%% Rerun at earliest convenience
corrAll(1,:) = mean(outputReconstructionMultipleROIs_E22.corrAll,1);
corrAll(2,:) = mean(outputReconstructionMultipleROIs_E39.corrAll,1);
corrAll(3,:) = mean(outputReconstructionMultipleROIs_E43.corrAll,1);
corrAll(4,:) = mean(outputReconstructionMultipleROIs_E44.corrAll,1);
corrAll(5,:) = mean(outputReconstructionMultipleROIs_E47.corrAll,1);
corrAll(6,:) = mean(outputReconstructionMultipleROIs_E48.corrAll,1);
corrAll(7,:) = mean(outputReconstructionMultipleROIs_E65.corrAll,1);

figure;
nieh_barSEM(corrAll);
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
