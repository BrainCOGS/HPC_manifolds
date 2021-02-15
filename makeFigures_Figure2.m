
%% Load the fnameStruct variable to get file names

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');


%% Make Choice Specific Sequence Plots

% Make sure to be in the shuffles folder, otherwise it will regenerate the
% shuffles, which will take some time

metamouse_Y         = generateMetamouse([5 2], 'noLog', 1, 'keepTrials', 100, ...
    fnameStruct, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, ...
    {[], []}, 'towers', 'all', 'both', 1);
 
plot_metamouse_seqplot(metamouse_Y);

 
 %% Make Evidence Sequence Plot
 
metamouse_E = generateMetamouse([5 2], 'noLog', 1, 'keepTrials', 100, ...
    fnameStruct, 2, {'Evidence'}, {[-15:16]}, {'Evidence', 'Position'}, ...
    {[], []}, 'towers', 'all', 'both', 0);
 
plot_metamouse_seqplot(metamouse_E);

% This also generates the training/testing plots in Extended Data Figure 1


%% Calculate Position x Evidence Skaggs for all animals
 
out_E22_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_ExY_[5_2]_%d.mat', 1, fnameStruct(1).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
out_E39_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_ExY_[5_2]_%d.mat', 1, fnameStruct(2).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
out_E43_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_ExY_[5_2]_%d.mat', 1, fnameStruct(3).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
out_E44_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_ExY_[5_2]_%d.mat', 1, fnameStruct(4).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
out_E47_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_ExY_[5_2]_%d.mat', 1, fnameStruct(5).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
out_E48_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_ExY_[5_2]_%d.mat', 1, fnameStruct(6).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');
out_E65_ExY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_ExY_[5_2]_%d.mat', 1, fnameStruct(7).fname, 2, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both');

%% Plot the ExY examples for E44

plot_getSkaggs_summarySLIM(out_E44_ExY, [325 8 11 28 30 37 38 76 77 79 97 104 116 123 160 176 183 231 247 296 1 378 435 465 617],25,[5 5]);

%% Run the random Dimension Data

rng(1);

out_E22_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
out_E39_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
out_E43_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
out_E44_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
out_E47_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
out_E48_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
out_E65_ExR = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');

out_E22_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E22_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
out_E39_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E39_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
out_E43_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E43_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
out_E44_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E44_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
out_E47_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E47_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
out_E48_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E48_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
out_E65_RxY = getSkaggs([5 2], 'noLog', 1, 'keepTrials', 100, 'E65_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');


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

figure; 
subplot(1,2,1)
bar(skaggsMetric_All_mean_Sig);
hold on;
er = errorbar([1:3], skaggsMetric_All_mean_Sig, [0 0 0], skaggsMetric_All_sem_Sig);
er.LineStyle = 'none'; 
er.Color = 'b';   
ylim([0 .03]);
xticklabels({'E x Y','Random x Y','E x Random'});
xtickangle(45)
ylabel('Mean Mutual Information Value');
title('Significant ROIs');

subplot(1,2,2)
bar(skaggsMetric_All_mean);
hold on;
er = errorbar([1:3], skaggsMetric_All_mean, [0 0 0], skaggsMetric_All_sem);
er.LineStyle = 'none'; 
er.Color = 'b';   
ylim([0 .03]);
xticklabels({'E x Y','Random x Y','E x Random'});
xtickangle(45)
title('All ROIs');

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

% Multiply values by 3 because there are three comparisons
[p1, r1] = ttest(skaggsMetric_All_ExY_Real_Sig, skaggsMetric_All_RxY_Real_Sig);
r1=r1*3

[p2, r2] = ttest(skaggsMetric_All_ExY_Real_Sig, skaggsMetric_All_ExR_Real_Sig);
r2=r2*3

[p3, r3] = ttest(skaggsMetric_All_ExR_Real_Sig, skaggsMetric_All_RxY_Real_Sig);
r3=r3*3