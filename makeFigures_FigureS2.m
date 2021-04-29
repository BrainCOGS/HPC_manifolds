
%% The individual field examples are generated with the code for the other
%% Related panels in Figure 2

%% Count ExR, RxY, and ExY

% First load the ExR, RxY, and ExY Data generated from Figure 2
load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExY_all.mat')
load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExR_and_RxY_all.mat')

Yonly(1)    = length(setdiff(out_E22_RxY.skaggsMetric.sigROIs, out_E22_ExY.skaggsMetric.sigROIs));
Eonly(1)    = length(setdiff(out_E22_ExR.skaggsMetric.sigROIs, out_E22_ExY.skaggsMetric.sigROIs));
ExYboth(1)  = length(out_E22_ExY.skaggsMetric.sigROIs);
totalROI(1) = length(out_E22_ExY.skaggsMetric.skaggs_real);
Yall(1)     = length(out_E22_RxY.skaggsMetric.sigROIs);
Eall(1)     = length(out_E22_ExR.skaggsMetric.sigROIs);

Yonly(2)    = length(setdiff(out_E39_RxY.skaggsMetric.sigROIs, out_E39_ExY.skaggsMetric.sigROIs));
Eonly(2)    = length(setdiff(out_E39_ExR.skaggsMetric.sigROIs, out_E39_ExY.skaggsMetric.sigROIs));
ExYboth(2)  = length(out_E39_ExY.skaggsMetric.sigROIs);
totalROI(2) = length(out_E39_ExY.skaggsMetric.skaggs_real);
Yall(2)     = length(out_E39_RxY.skaggsMetric.sigROIs);
Eall(2)     = length(out_E39_ExR.skaggsMetric.sigROIs);

Yonly(3)    = length(setdiff(out_E43_RxY.skaggsMetric.sigROIs, out_E43_ExY.skaggsMetric.sigROIs));
Eonly(3)    = length(setdiff(out_E43_ExR.skaggsMetric.sigROIs, out_E43_ExY.skaggsMetric.sigROIs));
ExYboth(3)  = length(out_E43_ExY.skaggsMetric.sigROIs);
totalROI(3) = length(out_E43_ExY.skaggsMetric.skaggs_real);
Yall(3)     = length(out_E43_RxY.skaggsMetric.sigROIs);
Eall(3)     = length(out_E43_ExR.skaggsMetric.sigROIs);

Yonly(4)    = length(setdiff(out_E44_RxY.skaggsMetric.sigROIs, out_E44_ExY.skaggsMetric.sigROIs));
Eonly(4)    = length(setdiff(out_E44_ExR.skaggsMetric.sigROIs, out_E44_ExY.skaggsMetric.sigROIs));
ExYboth(4)  = length(out_E44_ExY.skaggsMetric.sigROIs);
totalROI(4) = length(out_E44_ExY.skaggsMetric.skaggs_real);
Yall(4)     = length(out_E44_RxY.skaggsMetric.sigROIs);
Eall(4)     = length(out_E44_ExR.skaggsMetric.sigROIs);

Yonly(5)    = length(setdiff(out_E47_RxY.skaggsMetric.sigROIs, out_E47_ExY.skaggsMetric.sigROIs));
Eonly(5)    = length(setdiff(out_E47_ExR.skaggsMetric.sigROIs, out_E47_ExY.skaggsMetric.sigROIs));
ExYboth(5)  = length(out_E47_ExY.skaggsMetric.sigROIs);
totalROI(5) = length(out_E47_ExY.skaggsMetric.skaggs_real);
Yall(5)     = length(out_E47_RxY.skaggsMetric.sigROIs);
Eall(5)     = length(out_E47_ExR.skaggsMetric.sigROIs);

Yonly(6)    = length(setdiff(out_E48_RxY.skaggsMetric.sigROIs, out_E48_ExY.skaggsMetric.sigROIs));
Eonly(6)    = length(setdiff(out_E48_ExR.skaggsMetric.sigROIs, out_E48_ExY.skaggsMetric.sigROIs));
ExYboth(6)  = length(out_E48_ExY.skaggsMetric.sigROIs);
totalROI(6) = length(out_E48_ExY.skaggsMetric.skaggs_real);
Yall(6)     = length(out_E48_RxY.skaggsMetric.sigROIs);
Eall(6)     = length(out_E48_ExR.skaggsMetric.sigROIs);

Yonly(7)    = length(setdiff(out_E65_RxY.skaggsMetric.sigROIs, out_E65_ExY.skaggsMetric.sigROIs));
Eonly(7)    = length(setdiff(out_E65_ExR.skaggsMetric.sigROIs, out_E65_ExY.skaggsMetric.sigROIs));
ExYboth(7)  = length(out_E65_ExY.skaggsMetric.sigROIs);
totalROI(7) = length(out_E65_ExY.skaggsMetric.skaggs_real);
Yall(7)     = length(out_E65_RxY.skaggsMetric.sigROIs);
Eall(7)     = length(out_E65_ExR.skaggsMetric.sigROIs);

Per_Yonly = Yonly./totalROI;
Per_Eonly = Eonly./totalROI;
Per_Yall  = Yall./totalROI;
Per_Eall  = Eall./totalROI;
Per_ExYBoth = ExYboth./totalROI;
Per_ExYInd  = Per_Yall.*Per_Eall;

figure;
nieh_barSEMpaired(Per_ExYBoth, Per_ExYInd);
sourceData_S2f = [Per_ExYBoth' Per_ExYInd'];
signrank(Per_ExYBoth, Per_ExYInd)
xticklabels({'p(ExY)', 'p(RxY)*p(ExR)'});

% Venn Diagram
figure;
venn([sum(ExYboth)+sum(Yonly) sum(ExYboth)+sum(Eonly)],sum(ExYboth));
sourceData_S2d = [sum(Yonly) sum(Eonly) sum(ExYboth)];
sum_Yonly   = sum(Yonly)/sum(totalROI)
sum_Eonly   = sum(Eonly)/sum(totalROI)
sum_ExYboth = sum(ExYboth)/sum(totalROI)
sum_none    = 1-(sum_Yonly+sum_Eonly+sum_ExYboth)


%% Find distribution of ExR and RxY skaggs values

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

load("C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExY_all.mat");

load("C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExRtemp_and_RxYtemp_FigS2_cutdown.mat")
% Or run the following
% Make sure you're in the shuffle folder!!

% rng(42);
% for i=1:50
%     out_E22_ExRtemp{i} = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E22_ExR_[5_2]_%d.mat', 1, fnameStruct(1).fname, 2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E22_RxYtemp{i} = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E22_RxY_[5_2]_%d.mat', 1, fnameStruct(1).fname, 2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E22(i,:) = out_E22_ExRtemp{i}.skaggsMetric.skaggs_real;
%     skaggsDRxY_E22(i,:) = out_E22_RxYtemp{i}.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E22: ' num2str(i) ' of 50 shuffles done']);
% end
% 
% rng(42);
% for i=1:50
%     out_E39_ExRtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E39_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E39_RxYtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E39_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat', 2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E39(i,:) = out_E39_ExRtemp.skaggsMetric.skaggs_real;
%     skaggsDRxY_E39(i,:) = out_E39_RxYtemp.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E39: ' num2str(i) ' of 50 shuffles done']);
% end
% 
% rng(42);
% for i=1:50
%     out_E43_ExRtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E43_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E43_RxYtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E43_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E43(i,:) = out_E43_ExRtemp.skaggsMetric.skaggs_real;
%     skaggsDRxY_E43(i,:) = out_E43_RxYtemp.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E43: ' num2str(i) ' of 50 shuffles done']);
% end
% 
% rng(42);
% for i=1:50
%     out_E44_ExRtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E44_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E44_RxYtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E44_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E44(i,:) = out_E44_ExRtemp.skaggsMetric.skaggs_real;
%     skaggsDRxY_E44(i,:) = out_E44_RxYtemp.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E44: ' num2str(i) ' of 50 shuffles done']);
% end
% 
% rng(42);
% for i=1:50
%     out_E47_ExRtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E47_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E47_RxYtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E47_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E47(i,:) = out_E47_ExRtemp.skaggsMetric.skaggs_real;
%     skaggsDRxY_E47(i,:) = out_E47_RxYtemp.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E47: ' num2str(i) ' of 50 shuffles done']);
% end
% 
% rng(42);
% for i=1:50
%     out_E48_ExRtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E48_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E48_RxYtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E48_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E48(i,:) = out_E48_ExRtemp.skaggsMetric.skaggs_real;
%     skaggsDRxY_E48(i,:) = out_E48_RxYtemp.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E48: ' num2str(i) ' of 50 shuffles done']);
% end
% 
% rng(42);
% for i=1:50
%     out_E65_ExRtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E65_ExR_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Evidence', 'Random'}, {[-15:16], []}, {'Position', 'Evidence'}, {[0:10:300], [-15:16]}, 'towers', 'all', 'both');
%     out_E65_RxYtemp = getSkaggs_cutdown([5 2], 'noLog', 1, 'keepTrials', 0, 'E65_RxY_[5_2]_%d.mat', 1, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  2, {'Position', 'Random'}, {[0:10:300], []}, {'Evidence', 'Position'}, {[-15:16], [0:10:300]}, 'towers', 'all', 'both');
%     skaggsDExR_E65(i,:) = out_E65_ExRtemp.skaggsMetric.skaggs_real;
%     skaggsDRxY_E65(i,:) = out_E65_RxYtemp.skaggsMetric.skaggs_real;
%     close all;
%     disp(['E65: ' num2str(i) ' of 50 shuffles done']);
% end

% Run the analysis
for i=1:7
    
    if i==1
        skaggsDExR = skaggsDExR_E22;
        skaggsDRxY = skaggsDRxY_E22;
        out_ExY    = out_E22_ExY;
    elseif i==2
        skaggsDExR = skaggsDExR_E39;
        skaggsDRxY = skaggsDRxY_E39;
        out_ExY    = out_E39_ExY;
    elseif i==3
        skaggsDExR = skaggsDExR_E43;
        skaggsDRxY = skaggsDRxY_E43;
        out_ExY    = out_E43_ExY;
    elseif i==4
        skaggsDExR = skaggsDExR_E44;
        skaggsDRxY = skaggsDRxY_E44;
        out_ExY    = out_E44_ExY;
    elseif i==5
        skaggsDExR = skaggsDExR_E47;
        skaggsDRxY = skaggsDRxY_E47;
        out_ExY    = out_E47_ExY;
    elseif i==6
        skaggsDExR = skaggsDExR_E48;
        skaggsDRxY = skaggsDRxY_E48;
        out_ExY    = out_E48_ExY;
    elseif i==7
        skaggsDExR = skaggsDExR_E65;
        skaggsDRxY = skaggsDRxY_E65;
        out_ExY    = out_E65_ExY;
    end
    
    skaggsDExR_mean = mean(skaggsDExR);
    skaggsDExR_std = std(skaggsDExR);
    skaggsDExR_thres = skaggsDExR_mean+(2*skaggsDExR_std);
    skaggsDExR_sig = out_ExY.skaggsMetric.skaggs_real>skaggsDExR_thres;
    skaggsROIAll(i).sigExR = sum(skaggsDExR_sig(out_ExY.skaggsMetric.sigROIs));
    sigROIDExR = find(skaggsDExR_sig(out_ExY.skaggsMetric.sigROIs)~=1);
    skaggsROIAll(i).sigROIDExR = sigROIDExR;
    
    skaggsDRxY_mean = mean(skaggsDRxY);
    skaggsDRxY_std = std(skaggsDRxY);
    skaggsDRxY_thres = skaggsDRxY_mean+(2*skaggsDRxY_std);
    skaggsDRxY_sig = out_ExY.skaggsMetric.skaggs_real>skaggsDRxY_thres;
    skaggsROIAll(i).sigRxY = sum(skaggsDRxY_sig(out_ExY.skaggsMetric.sigROIs));
    sigROIDRxY = find(skaggsDRxY_sig(out_ExY.skaggsMetric.sigROIs)~=1);
    skaggsROIAll(i).sigROIDRxY = sigROIDRxY;
    
    skaggsROIAll(i).lengthSig = length(out_ExY.skaggsMetric.sigROIs);
    
    skaggsROIAll(i).lengthROInotExR  = length([sigROIDExR]);
    skaggsROIAll(i).lengthROInotRxY  = length([sigROIDRxY]);
    
    skaggsROIAll(i).uniqueROInotBoth = unique([sigROIDExR sigROIDRxY]);
    skaggsROIAll(i).lengthROInotBoth = length(unique([sigROIDExR sigROIDRxY]));
    
end

sumGood  = sum([skaggsROIAll.lengthSig]) - sum([skaggsROIAll.lengthROInotBoth]);
sumFakeE = sum([skaggsROIAll.lengthROInotExR]);
sumFakeY = sum([skaggsROIAll.lengthROInotRxY]);
sumSig   = sum([skaggsROIAll.lengthSig]);

perGood = sumGood/sumSig*100
perFakeE = sumFakeE/sumSig*100
perFakeY = sumFakeY/sumSig*100

% Check that this is 1, nothing funky happening
sumGood+sumFakeE+sumFakeY == sum([skaggsROIAll.lengthSig])

figure;
labels = {'Good' ,'Fake E', 'Fake Y'};
pie([sumGood sumFakeE sumFakeY], labels);
sourceData_S2e = [sumFakeE sumFakeY sumGood sumSig];


%% Calculate number of fields

% First load the ExY Data generated from Figure 2
load("C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_ExY_all.mat");

pixelwiseAll_sig = [out_E22_ExY.pixelwise(out_E22_ExY.skaggsMetric.sigROIs) out_E39_ExY.pixelwise(out_E39_ExY.skaggsMetric.sigROIs) out_E43_ExY.pixelwise(out_E43_ExY.skaggsMetric.sigROIs) out_E44_ExY.pixelwise(out_E44_ExY.skaggsMetric.sigROIs) out_E47_ExY.pixelwise(out_E47_ExY.skaggsMetric.sigROIs) out_E48_ExY.pixelwise(out_E48_ExY.skaggsMetric.sigROIs) out_E65_ExY.pixelwise(out_E65_ExY.skaggsMetric.sigROIs)];
num_peaks_sig = peak_counterSLIM(pixelwiseAll_sig,9);
mean(num_peaks_sig)
nieh_sem(num_peaks_sig)

figure
histogram(num_peaks_sig, [-0.5:1:6.5], 'Normalization','probability')
sourceData_S2g = num_peaks_sig';
xlabel('number of peaks')
ylabel('frequency')
set(gca,'box','off')


