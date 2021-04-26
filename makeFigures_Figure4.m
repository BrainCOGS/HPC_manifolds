
%% Generate the doublets and triplets

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletsTriplets_all.mat')
% Or run this code

% outputDoublets_E22 = findDoublets(11,fnameStruct(1).fname,[5 0], 1);
% outputDoublets_E39 = findDoublets(11,fnameStruct(2).fname,[5 0], 1);
% outputDoublets_E43 = findDoublets(5, fnameStruct(3).fname,[5 0], 1);
% outputDoublets_E44 = findDoublets(5, fnameStruct(4).fname,[5 0], 1);
% outputDoublets_E47 = findDoublets(5, fnameStruct(5).fname,[5 0], 1);
% outputDoublets_E48 = findDoublets(5, fnameStruct(6).fname,[5 0], 1);
% outputDoublets_E65 = findDoublets(5, fnameStruct(7).fname,[5 0], 1);
% 
% outputTriplets_E22 = findTriplets(outputDoublets_E22);
% outputTriplets_E39 = findTriplets(outputDoublets_E39);
% outputTriplets_E43 = findTriplets(outputDoublets_E43);
% outputTriplets_E44 = findTriplets(outputDoublets_E44);
% outputTriplets_E47 = findTriplets(outputDoublets_E47);
% outputTriplets_E48 = findTriplets(outputDoublets_E48);
% outputTriplets_E65 = findTriplets(outputDoublets_E65);

%% Plot the example doublets (Fig. 4a)

load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletsTriplets_all.mat')
load(fnameStruct(2).fname);
deltaT = score.deltaT;

% ****************IMPORTANT***************************************
% Need to manually delete the zeros at the end because of how matlab
% generates the cell arrays
sourceData_4a_ex1 = plotSequenceTrialROI(219,outputDoublets_E39, fnameStruct(2).fname, deltaT);
sourceData_4a_ex2 = plotSequenceTrialROI(122,outputDoublets_E39, fnameStruct(2).fname, deltaT);


%% Make the metamouse and generate plots (Fig. 4b-c, f-j 

makeOutputDoublets_Metamouse
makeOutputTriplets_Metamouse

plotDoublets(outputDoublets_metamouse);
plotTriplets(outputTriplets_metamouse);


%% To plot the trajectory between two cells in doublets (Fig. 4d)

% takes about 1 min to run
outputTrajectory_E39 = mind_plotDoubletTrajectory(fnameStruct(2).fname_mani, fnameStruct(2).fname, 1, 11, [5 0], 3, outputDoublets_E39, 219, [6 7], 0);


%% Make scatter plot of distance between doublet events (Fig. 4e)

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Either regenerate doublets or load

load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletManiCorrs_all.mat')

% Or run this code:

% load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletsTriplets_all.mat')
% outputDoubletManiCorrs_E22 = calcDoubletManifoldCorrelation(fnameStruct(1).fname, fnameStruct(1).fname_mani, 100, 5, outputDoublets_E22);
% outputDoubletManiCorrs_E39 = calcDoubletManifoldCorrelation(fnameStruct(2).fname, fnameStruct(2).fname_mani, 100, 5, outputDoublets_E39);
% outputDoubletManiCorrs_E43 = calcDoubletManifoldCorrelation(fnameStruct(3).fname, fnameStruct(3).fname_mani, 100, 5, outputDoublets_E43);
% outputDoubletManiCorrs_E44 = calcDoubletManifoldCorrelation(fnameStruct(4).fname, fnameStruct(4).fname_mani, 100, 5, outputDoublets_E44);
% outputDoubletManiCorrs_E47 = calcDoubletManifoldCorrelation(fnameStruct(5).fname, fnameStruct(5).fname_mani, 100, 5, outputDoublets_E47);
% outputDoubletManiCorrs_E48 = calcDoubletManifoldCorrelation(fnameStruct(6).fname, fnameStruct(6).fname_mani, 100, 5, outputDoublets_E48);
% outputDoubletManiCorrs_E65 = calcDoubletManifoldCorrelation(fnameStruct(7).fname, fnameStruct(7).fname_mani, 100, 5, outputDoublets_E65);

for i=1:length(fnameStruct)
    load(fnameStruct(i).fname);
    deltaTAll(i) = score.deltaT;
end

% Not really necessary since deltaTAll is identical
deltaT = mean(deltaTAll);

meanCorrs(1,1) = nanmean(outputDoubletManiCorrs_E22.real_corrs);
meanCorrs(1,2) = mean(nanmean(outputDoubletManiCorrs_E22.shuffled_corrs,1));
meanCorrs(2,1) = nanmean(outputDoubletManiCorrs_E39.real_corrs);
meanCorrs(2,2) = mean(nanmean(outputDoubletManiCorrs_E39.shuffled_corrs,1));
meanCorrs(3,1) = nanmean(outputDoubletManiCorrs_E43.real_corrs);
meanCorrs(3,2) = mean(nanmean(outputDoubletManiCorrs_E43.shuffled_corrs,1));
meanCorrs(4,1) = nanmean(outputDoubletManiCorrs_E44.real_corrs);
meanCorrs(4,2) = mean(nanmean(outputDoubletManiCorrs_E44.shuffled_corrs,1));
meanCorrs(5,1) = nanmean(outputDoubletManiCorrs_E47.real_corrs);
meanCorrs(5,2) = mean(nanmean(outputDoubletManiCorrs_E47.shuffled_corrs,1));
meanCorrs(6,1) = nanmean(outputDoubletManiCorrs_E48.real_corrs);
meanCorrs(6,2) = mean(nanmean(outputDoubletManiCorrs_E48.shuffled_corrs,1));
meanCorrs(7,1) = nanmean(outputDoubletManiCorrs_E65.real_corrs);
meanCorrs(7,2) = mean(nanmean(outputDoubletManiCorrs_E65.shuffled_corrs,1));

signrank(meanCorrs(:,1), meanCorrs(:,2))

plotCorrs(:,1) = [outputDoubletManiCorrs_E22.dist_mani_all_vertcat; outputDoubletManiCorrs_E39.dist_mani_all_vertcat; outputDoubletManiCorrs_E43.dist_mani_all_vertcat; outputDoubletManiCorrs_E44.dist_mani_all_vertcat; outputDoubletManiCorrs_E47.dist_mani_all_vertcat; outputDoubletManiCorrs_E48.dist_mani_all_vertcat; outputDoubletManiCorrs_E65.dist_mani_all_vertcat];
plotCorrs(:,2) = [outputDoubletManiCorrs_E22.dist_time_all_vertcat; outputDoubletManiCorrs_E39.dist_time_all_vertcat; outputDoubletManiCorrs_E43.dist_time_all_vertcat; outputDoubletManiCorrs_E44.dist_time_all_vertcat; outputDoubletManiCorrs_E47.dist_time_all_vertcat; outputDoubletManiCorrs_E48.dist_time_all_vertcat; outputDoubletManiCorrs_E65.dist_time_all_vertcat];

figure; 
scatter(outputDoubletManiCorrs_E22.dist_mani_all_vertcat,outputDoubletManiCorrs_E22.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
hold on;
scatter(outputDoubletManiCorrs_E39.dist_mani_all_vertcat,outputDoubletManiCorrs_E39.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
scatter(outputDoubletManiCorrs_E43.dist_mani_all_vertcat,outputDoubletManiCorrs_E43.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
scatter(outputDoubletManiCorrs_E44.dist_mani_all_vertcat,outputDoubletManiCorrs_E44.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
scatter(outputDoubletManiCorrs_E47.dist_mani_all_vertcat,outputDoubletManiCorrs_E47.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
scatter(outputDoubletManiCorrs_E48.dist_mani_all_vertcat,outputDoubletManiCorrs_E48.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
scatter(outputDoubletManiCorrs_E65.dist_mani_all_vertcat,outputDoubletManiCorrs_E65.dist_time_all_vertcat, 5,'o','filled','MarkerFaceAlpha',.1)
xlabel('Distance in Manifold Space');
ylabel('Time (s)');
legend('E22','E39','E43','E44','E47','E48','E65');
axis square


%% Number of trials in which doublet is active

load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletsTriplets_all.mat')

perTrials_E22 = [outputDoublets_E22.saveAll_doublets.number_doublets]/length(outputDoublets_E22.saveAll_trials);
perTrials_E39 = [outputDoublets_E39.saveAll_doublets.number_doublets]/length(outputDoublets_E39.saveAll_trials);
perTrials_E43 = [outputDoublets_E43.saveAll_doublets.number_doublets]/length(outputDoublets_E43.saveAll_trials);
perTrials_E44 = [outputDoublets_E44.saveAll_doublets.number_doublets]/length(outputDoublets_E44.saveAll_trials);
perTrials_E47 = [outputDoublets_E47.saveAll_doublets.number_doublets]/length(outputDoublets_E47.saveAll_trials);
perTrials_E48 = [outputDoublets_E48.saveAll_doublets.number_doublets]/length(outputDoublets_E48.saveAll_trials);
perTrials_E65 = [outputDoublets_E65.saveAll_doublets.number_doublets]/length(outputDoublets_E65.saveAll_trials);

perTrials_ALL = [perTrials_E22 perTrials_E39 perTrials_E43 perTrials_E44 perTrials_E47 perTrials_E48 perTrials_E65];
mean(perTrials_ALL)*100
nieh_sem(perTrials_ALL)*100


%% ROC Analysis

load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletsTriplets_all.mat')

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

load('C:\Neuroscience\imaging\FINAL\roc_Data\outputOptimalROC_allAnimals.mat')

% Or run this code:

% outputOptimalROC_E22_5 = mind_calcOptimalROC(fnameStruct(1).fname, fnameStruct(1).fname_mani, outputDoublets_E22,5);
% outputOptimalROC_E39_5 = mind_calcOptimalROC(fnameStruct(2).fname, fnameStruct(2).fname_mani, outputDoublets_E39,5);
% outputOptimalROC_E43_5 = mind_calcOptimalROC(fnameStruct(3).fname, fnameStruct(3).fname_mani, outputDoublets_E43,5);
% outputOptimalROC_E44_5 = mind_calcOptimalROC(fnameStruct(4).fname, fnameStruct(4).fname_mani, outputDoublets_E44,5);
% outputOptimalROC_E47_5 = mind_calcOptimalROC(fnameStruct(5).fname, fnameStruct(5).fname_mani, outputDoublets_E47,5);
% outputOptimalROC_E48_5 = mind_calcOptimalROC(fnameStruct(6).fname, fnameStruct(6).fname_mani, outputDoublets_E48,5);
% outputOptimalROC_E65_5 = mind_calcOptimalROC(fnameStruct(7).fname, fnameStruct(7).fname_mani, outputDoublets_E65,5);

ValuesROC_ALL(1,:) = outputOptimalROC_E22_5.ValuesROC;
ValuesROC_ALL(2,:) = outputOptimalROC_E39_5.ValuesROC;
ValuesROC_ALL(3,:) = outputOptimalROC_E43_5.ValuesROC;
ValuesROC_ALL(4,:) = outputOptimalROC_E44_5.ValuesROC;
ValuesROC_ALL(5,:) = outputOptimalROC_E47_5.ValuesROC;
ValuesROC_ALL(6,:) = outputOptimalROC_E48_5.ValuesROC;
ValuesROC_ALL(7,:) = outputOptimalROC_E65_5.ValuesROC;

mean(ValuesROC_ALL)
nieh_sem(ValuesROC_ALL)

figure;
hold on;
plot(outputOptimalROC_E22_5.totalROC_ALL(:,2), outputOptimalROC_E22_5.totalROC_ALL(:,1), 'LineWidth',1);
plot(outputOptimalROC_E39_5.totalROC_ALL(:,2), outputOptimalROC_E39_5.totalROC_ALL(:,1), 'LineWidth',1);
plot(outputOptimalROC_E43_5.totalROC_ALL(:,2), outputOptimalROC_E43_5.totalROC_ALL(:,1), 'LineWidth',1);
plot(outputOptimalROC_E44_5.totalROC_ALL(:,2), outputOptimalROC_E44_5.totalROC_ALL(:,1), 'LineWidth',1);
plot(outputOptimalROC_E47_5.totalROC_ALL(:,2), outputOptimalROC_E47_5.totalROC_ALL(:,1), 'LineWidth',1);
plot(outputOptimalROC_E48_5.totalROC_ALL(:,2), outputOptimalROC_E48_5.totalROC_ALL(:,1), 'LineWidth',1);
plot(outputOptimalROC_E65_5.totalROC_ALL(:,2), outputOptimalROC_E65_5.totalROC_ALL(:,1), 'LineWidth',1);
axis square
xlabel('False Positive Rate');
ylabel('True Positive Rate');
xlim([0 1]);
ylim([0 1]);
xlim([0 1]);
xticks([0:.2:1]);
ylim([0 1]);
yticks([0:.2:1]);
xlabel('False Positive Rate');
ylabel('True Positive Rate');

legend('E22', 'E39', 'E43', 'E44', 'E47', 'E48', 'E65');


