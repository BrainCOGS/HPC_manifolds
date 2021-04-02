
%% Intrinsic dimensionality

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

outputFitPowerLaw = mind_fitPowerLaw_FunctionSLIM(fnameStruct, 'towers', 1, 5);

% Saved the output as outputFitPowerLaw_Towers to use with the T7 data in S3
% Location: "C:\Neuroscience\imaging\FINAL\fitPowerLaw_Data\outputFitPowerLaw_Towers.mat"

%% Example reconstruction

% Copy and paste this to matlab on spock to submit all the jobs
% 150 is not there because it is a bad trial

% trialsList = [145 146 147 148 149 151 152 153 154 157 158 159 160];
% 
% for i=1:length(trialsList)
%     curTrial = num2str(trialsList(i));
%     command_string = redo_deadjob_function(curTrial);
%     submitJobs(command_string, 'towers');
% end

%after jobs have completed, run:
collect_allTrials_reconstructionMultipleTrials


%% Reconstruction scores

clear all; 

% First run run_analysis_allTrials on matlab on spock to generate the
% held-out trial data

collect_allTrials_2to7_minLeaves500

figure; 
nieh_barSEM(maxReconstruct(1,:), maxReconstruct(2,:), maxReconstruct(3,:), maxReconstruct(4,:), maxReconstruct(5,:), maxReconstruct(6,:));
ylabel('Mean Cross-Validated Score');
xticklabels({'2', '3', '4', '5', '6', '7'});
xlabel('Num Embedding Dims');
set(gca, 'box', 'off')
axis square

%% Comparison with PCA

% Need the data from collect_allTrials_2to7_minLeaves500 from above

mean_trials_animals = mind_pcaTest(fnameStruct);

meanMaxReconstruct = mean(maxReconstruct,2);

% Then find the index that corresponds to the first entry that surpasses
% the numbers in the reconstruction scores plot above.

greater4 = mean_trials_animals>meanMaxReconstruct(3);
greater4 = find(greater4==1);
greater4 = greater4(1);
greater5 = mean_trials_animals>meanMaxReconstruct(4);
greater5 = find(greater5==1);
greater5 = greater5(1);
greater6 = mean_trials_animals>meanMaxReconstruct(5);
greater6 = find(greater6==1);
greater6 = greater6(1);

[greater4 greater5 greater6]


%% Tiled fields

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Important functions it uses are:
% mind_fitFiringFieldsNEW_dimX_manuel, extractVariables
outputTiledFields = mind_plotTiledFieldsSLIM(fnameStruct(7).fname, fnameStruct(7).fname_mani);

% To generate the movie, after running mind_plotTiledFields:
movieName = 'Supp_movie_2';
mind_makeTiledMovie(outputTiledFields.manifold3d, outputTiledFields.ROIactivities_thres, movieName);


%% Position and evidence gradients

load(fnameStruct(7).fname_mani);
mind_plotManifoldGradients(outMind, fnameStruct(7).fname, 'towers',1)
set(gcf,'renderer','painters');


%% Decode position and evidence

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Load data from "C:\Neuroscience\imaging\FINAL\decoding_Data\decodeEandY_all.mat"
% Or run the following

dimEmbedList = [2:7];
for i=1:length(fnameStruct)
    for j=1:length(dimEmbedList)
        outputNonlinearDecoding_E = mind_nonlinearDecoding_dimX_All(fnameStruct(i).fname, fnameStruct(i).fname_mani,5,'GP','Evidence','towers',0,1,[], dimEmbedList(j));
        outputNonlinearDecodingAll_E{i,j} = outputNonlinearDecoding_E;
        meancorrAll_E(i,j) = outputNonlinearDecoding_E.meancorr;
        
        outputNonlinearDecoding_Y = mind_nonlinearDecoding_dimX_All(fnameStruct(i).fname, fnameStruct(i).fname_mani,5,'GP','Position','towers',0,1,[], dimEmbedList(j));
        outputNonlinearDecodingAll_Y{i,j} = outputNonlinearDecoding_Y;
        meancorrAll_Y(i,j) = outputNonlinearDecoding_Y.meancorr;
        
        disp(['Animal ' num2str(i) ' of ' num2str(length(fnameStruct)) ', dim ' num2str(dimEmbedList(j)) ' finished']);
    end
    
    outputNonlinearDecoding_ROIs_E = mind_nonlinearDecoding_dimX_All(fnameStruct(i).fname, fnameStruct(i).fname_mani,5,'GP','Evidence','towers',0,0,[],[]);
    outputNonlinearDecoding_ROIsAll_E{i,j} = outputNonlinearDecoding_ROIs_E;
    meancorrROIsAll_E(i) = outputNonlinearDecoding_ROIs_E.meancorr;
    
    outputNonlinearDecoding_ROIs_Y = mind_nonlinearDecoding_dimX_All(fnameStruct(i).fname, fnameStruct(i).fname_mani,5,'GP','Position','towers',0,0,[],[]);
    outputNonlinearDecoding_ROIsAll_Y{i,j} = outputNonlinearDecoding_ROIs_Y;
    meancorrROIsAll_Y(i) = outputNonlinearDecoding_ROIs_Y.meancorr;
    
    disp(['Animal ' num2str(i) ' of ' num2str(length(fnameStruct)) ', ROI finished']);
end

figure; 
subplot(1,2,1)
nieh_barSEM(meancorrAll_Y, meancorrROIsAll_Y);
hold on;
        scatter([ones(length(fnameStruct),1); ...
                 ones(length(fnameStruct),1)*2; ...
                 ones(length(fnameStruct),1)*3; ...
                 ones(length(fnameStruct),1)*4; ...
                 ones(length(fnameStruct),1)*5; ...
                 ones(length(fnameStruct),1)*6; ...
                 ones(length(fnameStruct),1)*7; ...
                 ],[meancorrAll_Y(:); meancorrROIsAll_Y']);
ylabel('Decoding Index (r)');
xticklabels({'2', '3', '4', '5', '6', '7', 'ROIs'});
xlabel('# Dims embedded, last bar is ROIs');
title('Position')
set(gca, 'box', 'off')

subplot(1,2,2)
nieh_barSEM(meancorrAll_E, meancorrROIsAll_E);
hold on;
        scatter([ones(length(fnameStruct),1); ...
                 ones(length(fnameStruct),1)*2; ...
                 ones(length(fnameStruct),1)*3; ...
                 ones(length(fnameStruct),1)*4; ...
                 ones(length(fnameStruct),1)*5; ...
                 ones(length(fnameStruct),1)*6; ...
                 ones(length(fnameStruct),1)*7; ...
                 ],[meancorrAll_E(:); meancorrROIsAll_E']);
ylabel('Decoding Index (r)');
xticklabels({'2', '3', '4', '5', '6', '7', 'ROIs'});
xlabel('# Dims embedded, last bar is ROIs');
title('Evidence')
set(gca, 'box', 'off')


%% Align Multiple Animals (Fig. 3i, j)

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Run the code above or
load("C:\Neuroscience\imaging\FINAL\decoding_Data\decodeEandY_all.mat")

% index of 4 means 5 dim
bestfit_position = meancorrAll_Y(:,4);
bestfit_evidence = meancorrAll_E(:,4);

load('C:\Neuroscience\imaging\FINAL\HPC2HPC_Data\outputHPC2HPC.mat')
% Or run this code:
%outputHPC2HPC = map_HPC2HPC(fnameStruct, bestfit_position, bestfit_evidence);

outputPlotHPC2HPC = plot_HPC2HPC(outputHPC2HPC);

