
%% Load the fnameStruct variable to get file names

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');


%% Calculate Position and Evidence skaggs for all animals (used later)

load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_Y_all.mat')
load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_E_all.mat')

% or run this code:
% out_E22_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_Y_[5_2]_%d.mat', 1, fnameStruct(1).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E39_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_Y_[5_2]_%d.mat', 1, fnameStruct(2).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E43_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_Y_[5_2]_%d.mat', 1, fnameStruct(3).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E44_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_Y_[5_2]_%d.mat', 1, fnameStruct(4).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E47_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_Y_[5_2]_%d.mat', 1, fnameStruct(5).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E48_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_Y_[5_2]_%d.mat', 1, fnameStruct(6).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E65_Y = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_Y_[5_2]_%d.mat', 1, fnameStruct(7).fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% 
% out_E22_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_E_[5_2]_%d.mat', 1, fnameStruct(1).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E39_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_E_[5_2]_%d.mat', 1, fnameStruct(2).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E43_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_E_[5_2]_%d.mat', 1, fnameStruct(3).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E44_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_E_[5_2]_%d.mat', 1, fnameStruct(4).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E47_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_E_[5_2]_%d.mat', 1, fnameStruct(5).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E48_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_E_[5_2]_%d.mat', 1, fnameStruct(6).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
% out_E65_E = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_E_[5_2]_%d.mat', 1, fnameStruct(7).fname, 2, {'Evidence'}, {[-15:16]},   {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');


%% Make Choice Specific Sequence Plots

% Make sure to be in the shuffles folder, otherwise it will regenerate the
% shuffles, which will take some time

load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\metamouse_Y.mat')

% or run this code
% metamouse_Y         = generateMetamouse([5 2], 'noLog', 1, 'keepTrials', 100, ...
%     fnameStruct, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, ...
%     {[], []}, 'towers', 'all', 'both', 1);
 
plot_metamouse_seqplot(metamouse_Y);

 
 %% Make Evidence Sequence Plot
 
load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\metamouse_E.mat')

% or run this code
% metamouse_E = generateMetamouse([5 2], 'noLog', 1, 'keepTrials', 100, ...
%     fnameStruct, 2, {'Evidence'}, {[-15:16]}, {'Evidence', 'Position'}, ...
%     {[], []}, 'towers', 'all', 'both', 0);
 
plot_metamouse_seqplot(metamouse_E);

% This also generates the training/testing plots in Extended Data Figure 1


%% Calculate Position x Evidence Skaggs for all animals
 
load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExY_all.mat')

% or run this code
% out_E22_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_ExY_[5_2]_%d.mat', 1, fnameStruct(1).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
% out_E39_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_ExY_[5_2]_%d.mat', 1, fnameStruct(2).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
% out_E43_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_ExY_[5_2]_%d.mat', 1, fnameStruct(3).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
% out_E44_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_ExY_[5_2]_%d.mat', 1, fnameStruct(4).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
% out_E47_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_ExY_[5_2]_%d.mat', 1, fnameStruct(5).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
% out_E48_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_ExY_[5_2]_%d.mat', 1, fnameStruct(6).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
% out_E65_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_ExY_[5_2]_%d.mat', 1, fnameStruct(7).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');

%% Plot the ExY examples for E44

plot_getSkaggs_summarySLIM(out_E44_ExY, [325 8 11 28 30 37 38 76 77 79 97 104 116 123 160 176 183 231 247 296 1 378 435 465 617],25,[5 5]);

%% Run the random Dimension Data

load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExR_and_RxY_all.mat')

% or run this code

% rng(1);
% 
% out_E22_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% out_E39_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% out_E43_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% out_E44_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% out_E47_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% out_E48_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% out_E65_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
% 
% out_E22_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
% out_E39_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
% out_E43_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
% out_E44_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
% out_E47_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
% out_E48_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
% out_E65_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');


%% Calculate whether distribution of skagg values is more for real or shuffled dimension

skaggsMetric_All_ExY_Real_Sig = [out_E22_ExY.skaggsMetric.skaggs_real(out_E22_ExY.skaggsMetric.sigROIs) out_E39_ExY.skaggsMetric.skaggs_real(out_E39_ExY.skaggsMetric.sigROIs) out_E43_ExY.skaggsMetric.skaggs_real(out_E43_ExY.skaggsMetric.sigROIs) out_E44_ExY.skaggsMetric.skaggs_real(out_E44_ExY.skaggsMetric.sigROIs) out_E47_ExY.skaggsMetric.skaggs_real(out_E47_ExY.skaggsMetric.sigROIs) out_E48_ExY.skaggsMetric.skaggs_real(out_E48_ExY.skaggsMetric.sigROIs) out_E65_ExY.skaggsMetric.skaggs_real(out_E65_ExY.skaggsMetric.sigROIs)];
skaggsMetric_All_RxY_Real_Sig = [out_E22_RxY.skaggsMetric.skaggs_real(out_E22_ExY.skaggsMetric.sigROIs) out_E39_RxY.skaggsMetric.skaggs_real(out_E39_ExY.skaggsMetric.sigROIs) out_E43_RxY.skaggsMetric.skaggs_real(out_E43_ExY.skaggsMetric.sigROIs) out_E44_RxY.skaggsMetric.skaggs_real(out_E44_ExY.skaggsMetric.sigROIs) out_E47_RxY.skaggsMetric.skaggs_real(out_E47_ExY.skaggsMetric.sigROIs) out_E48_RxY.skaggsMetric.skaggs_real(out_E48_ExY.skaggsMetric.sigROIs) out_E65_RxY.skaggsMetric.skaggs_real(out_E65_ExY.skaggsMetric.sigROIs)];
skaggsMetric_All_ExR_Real_Sig = [out_E22_ExR.skaggsMetric.skaggs_real(out_E22_ExY.skaggsMetric.sigROIs) out_E39_ExR.skaggsMetric.skaggs_real(out_E39_ExY.skaggsMetric.sigROIs) out_E43_ExR.skaggsMetric.skaggs_real(out_E43_ExY.skaggsMetric.sigROIs) out_E44_ExR.skaggsMetric.skaggs_real(out_E44_ExY.skaggsMetric.sigROIs) out_E47_ExR.skaggsMetric.skaggs_real(out_E47_ExY.skaggsMetric.sigROIs) out_E48_ExR.skaggsMetric.skaggs_real(out_E48_ExY.skaggsMetric.sigROIs) out_E65_ExR.skaggsMetric.skaggs_real(out_E65_ExY.skaggsMetric.sigROIs)];

skaggsMetric_All_ExY_Real = [out_E22_ExY.skaggsMetric.skaggs_real out_E39_ExY.skaggsMetric.skaggs_real out_E43_ExY.skaggsMetric.skaggs_real out_E44_ExY.skaggsMetric.skaggs_real out_E47_ExY.skaggsMetric.skaggs_real out_E48_ExY.skaggsMetric.skaggs_real out_E65_ExY.skaggsMetric.skaggs_real];
skaggsMetric_All_RxY_Real = [out_E22_RxY.skaggsMetric.skaggs_real out_E39_RxY.skaggsMetric.skaggs_real out_E43_RxY.skaggsMetric.skaggs_real out_E44_RxY.skaggsMetric.skaggs_real out_E47_RxY.skaggsMetric.skaggs_real out_E48_RxY.skaggsMetric.skaggs_real out_E65_RxY.skaggsMetric.skaggs_real];
skaggsMetric_All_ExR_Real = [out_E22_ExR.skaggsMetric.skaggs_real out_E39_ExR.skaggsMetric.skaggs_real out_E43_ExR.skaggsMetric.skaggs_real out_E44_ExR.skaggsMetric.skaggs_real out_E47_ExR.skaggsMetric.skaggs_real out_E48_ExR.skaggsMetric.skaggs_real out_E65_ExR.skaggsMetric.skaggs_real];

skaggsMetric_All_mean_Sig = [mean(skaggsMetric_All_ExY_Real_Sig) mean(skaggsMetric_All_RxY_Real_Sig) mean(skaggsMetric_All_ExR_Real_Sig)];
skaggsMetric_All_sem_Sig  = [nieh_sem(skaggsMetric_All_ExY_Real_Sig) nieh_sem(skaggsMetric_All_RxY_Real_Sig) nieh_sem(skaggsMetric_All_ExR_Real_Sig)];

skaggsMetric_All_mean = [mean(skaggsMetric_All_ExY_Real) mean(skaggsMetric_All_RxY_Real) mean(skaggsMetric_All_ExR_Real)];
skaggsMetric_All_sem  = [nieh_sem(skaggsMetric_All_ExY_Real) nieh_sem(skaggsMetric_All_RxY_Real) nieh_sem(skaggsMetric_All_ExR_Real)];

%pair1 = log([skaggsMetric_All_ExY_Real_Sig-skaggsMetric_All_RxY_Real_Sig]);
%pair2 = log([skaggsMetric_All_ExY_Real_Sig-skaggsMetric_All_ExR_Real_Sig]);
%pair3 = log([skaggsMetric_All_RxY_Real_Sig-skaggsMetric_All_ExR_Real_Sig]);

skaggsMetric_All_ExY_Real_Sig_log = log(skaggsMetric_All_ExY_Real_Sig);
skaggsMetric_All_RxY_Real_Sig_log = log(skaggsMetric_All_RxY_Real_Sig);
skaggsMetric_All_ExR_Real_Sig_log = log(skaggsMetric_All_ExR_Real_Sig);

sourceData_2e = [skaggsMetric_All_ExY_Real_Sig_log; skaggsMetric_All_RxY_Real_Sig_log; skaggsMetric_All_ExR_Real_Sig_log]';

figure; 
%scatter(ones(size(skaggsMetric_All_ExY_Real_Sig)).*(1+(rand(size(skaggsMetric_All_ExY_Real_Sig))-0.5)/5),skaggsMetric_All_ExY_Real_Sig,5,'k','filled','MarkerFaceAlpha',0.2);
%scatter(ones(size(skaggsMetric_All_RxY_Real_Sig)).*(2+(rand(size(skaggsMetric_All_RxY_Real_Sig))-0.5)/5),skaggsMetric_All_RxY_Real_Sig,5,'k','filled','MarkerFaceAlpha',0.2);
%scatter(ones(size(skaggsMetric_All_ExR_Real_Sig)).*(3+(rand(size(skaggsMetric_All_ExR_Real_Sig))-0.5)/5),skaggsMetric_All_ExR_Real_Sig,5,'k','filled','MarkerFaceAlpha',0.2);
plot([skaggsMetric_All_ExY_Real_Sig_log; skaggsMetric_All_RxY_Real_Sig_log; skaggsMetric_All_ExR_Real_Sig_log],'k','LineWidth',.01)
hold on;
h=boxplot([skaggsMetric_All_ExY_Real_Sig_log; skaggsMetric_All_RxY_Real_Sig_log; skaggsMetric_All_ExR_Real_Sig_log]','Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','b');
set(h,{'linew'},{2})
%ylim([0 .15]);
xticklabels({'E x Y','R x Y','E x R'});
xtickangle(45)
ylabel('Log Mutual Information');
title('Significant ROIs');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [100, 340, 313, 420])

% Multiply values by 3 because there are three comparisons
[p1, r1] = ttest(skaggsMetric_All_ExY_Real_Sig, skaggsMetric_All_RxY_Real_Sig);
r1=r1*3

[p2, r2] = ttest(skaggsMetric_All_ExY_Real_Sig, skaggsMetric_All_ExR_Real_Sig);
r2=r2*3

[p3, r3] = ttest(skaggsMetric_All_ExR_Real_Sig, skaggsMetric_All_RxY_Real_Sig);
r3=r3*3


%% This was moved to S2

sourceData_S2b = [skaggsMetric_All_ExY_Real_Sig; skaggsMetric_All_RxY_Real_Sig]';
sourceData_S2c = [skaggsMetric_All_ExY_Real_Sig; skaggsMetric_All_ExR_Real_Sig]';

figure; 
subplot(1,2,1)
scatter(skaggsMetric_All_ExY_Real_Sig, skaggsMetric_All_RxY_Real_Sig, 'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
xlim([0 .15])
ylim([0 .15])
axis square
hold on
plot([0 .15], [0 .15]);
xlabel('ExY Space');
ylabel('RxY');

subplot(1,2,2)
scatter(skaggsMetric_All_ExY_Real_Sig, skaggsMetric_All_ExR_Real_Sig, 'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3);
xlim([0 .15])
ylim([0 .15])
axis square
hold on
plot([0 .15], [0 .15]);
xlabel('ExY Space');
ylabel('ExR');

suptitle('SkaggsMetric, All - Real & Sig');


