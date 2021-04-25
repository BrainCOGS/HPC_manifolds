
%% Intrinsic dimensionality for T7, one-side cues task

fnameStruct_T7 = mind_makeFnameStruct('Edward','T7', 'none');

load("C:\Neuroscience\imaging\FINAL\fitPowerLaw_Data\outputFitPowerLaw_T7.mat");
% Or run this
% outputFitPowerLaw_T7 = mind_fitPowerLaw_FunctionSLIM(fnameStruct_T7, 'T7', 1, 4);

load("C:\Neuroscience\imaging\FINAL\fitPowerLaw_Data\outputFitPowerLaw_Towers.mat")

ranksumDim = ranksum(outputFitPowerLaw.exponents_l, outputFitPowerLaw_T7.exponents_l)

sourceData_S3b_1 = outputFitPowerLaw.exponents_l';
sourceData_S3b_2 = outputFitPowerLaw_T7.exponents_l';

figure; 
nieh_barSEM(outputFitPowerLaw.exponents_l, outputFitPowerLaw_T7.exponents_l);
ylabel('Dimension Estimate')
xticklabels({'Towers', 'T7'});
title(['ranksum is: ' num2str(ranksumDim)]);
hold on;
scatter(ones(1,7),outputFitPowerLaw.exponents_l,'.')
scatter(ones(1,4)*2,outputFitPowerLaw_T7.exponents_l,'.')
set(gca, 'box', 'off')


%% For the choice-specific sequences

load("C:\Neuroscience\imaging\FINAL\getSkaggs_Data\metamouse_Y_T7.mat");

% Or run this
% metamouse_Y_T7 = generateMetamouseT7([5 2], 'noLog', 1, 'keepTrials', 100, ...
%     fnameStruct_T7, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, ...
%     {[], []}, 'towers', 'all', 'both', 1);
 
plot_metamouse_seqplot(metamouse_Y_T7);

