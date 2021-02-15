function outputDoubletManiCorrs = calcDoubletManifoldCorrelation(fname, fname_mani, numShuf, numDim, outputDoublets)
%
%  This function calculates the correlation of the observed time between
%  doublets, and their distance on the manifold. It returns these values
%  together with a random control, that is constructed by taking the same
%  time difference, but for random sections of the trajectory on the
%  manifold.
%
%  Input:  fname:          path to the modeling file
%          fname_mani:     path to the manifold
%          outputDoublets: struct from Edwards doublet finder code
%          shuffles:       number of shuffles.
%          numDim:         number of dimentions of manifold
%  Output: shuffled_corrs: [N_doublets X shuffles]: The random control.
%          real_corrs:     [N_doublets X 1]: Observed correlation.
%
%  Example use:
%  fname = '../../data/E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat';
%  fname_mani = '../../data/E65_20180202_downsample1_RNG1_leaves200_2_7.mat';
%  load('../../data/E65_sequences_Easy.mat')
%  [shuffled_corrs, real_corrs] = calculate_doublet_mani_correlation(fname, fname_mani, 100, 3, outputDoublets_E65)
%

rng(42);

argins.numShuf = numShuf;
argins.numDim  = numDim;

nic_output = extractVariables('all', 2, 'keepTrials', 1, 0, 0, 4, [0 0], fname,'none','towers',1,1);
behavioralVariables = nic_output.behavioralVariables;
trialn = behavioralVariables.Trial;
[data_trial_nrs,~,ic] = unique(trialn);
ROIactivities = outputDoublets.saveAll_basics.ROIactivities;
time1 = behavioralVariables.Time;
load(fname_mani);

%% Get the manifold and the doublets in nice shape for further processing.
sample = outMind.dataDFF;
pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.dat.mindparameters.pca.n);
manifold3d = outMind.dat.allembed([outMind.dat.allembed.d]==numDim).f2m.map.transform(pca_coords);
padSpots = find(outMind.Datarange==0);
manifold3d_pad = zeros(length(outMind.Datarange),numDim);
manifold3d_pad(outMind.Datarange,:) = manifold3d;
for i=1:length(padSpots)
    manifold3d_pad(padSpots(i),:) = manifold3d_pad(padSpots(i)-1,:);
end
manifold3d = manifold3d_pad;
lengthMani = size(manifold3d,1);
saveAll_doublets   = outputDoublets.saveAll_doublets;

%% Loop though all doublets, find deltaT and deltaL on manifold

corrs_per_doublet = zeros(length(saveAll_doublets), 1);
corrs_per_doublet_rnd = zeros(length(saveAll_doublets), numShuf);

f = waitbar(0);
for i=1:length(saveAll_doublets)
    waitbar(i/length(saveAll_doublets));
    
    curDoublet = saveAll_doublets(i).cells;
    cell1 = curDoublet(1);
    cell2 = curDoublet(2);
    trialList = saveAll_doublets(i).trials_appear;
    
    % for this doublet, find all trials where it occured
    dist_time = zeros(length(trialList),1);
    dist_mani = zeros(length(trialList),1);
    dist_mani_rnd = zeros(length(trialList),numShuf);
    
    for j=1:length(trialList)
        curTrial = data_trial_nrs(trialList(j));
        cell1time = find(trialn==curTrial & ROIactivities(:,cell1)==1);
        cell2time = find(trialn==curTrial & ROIactivities(:,cell2)==1);
        t0 = min(cell1time) - min(find(trialn==curTrial));
        timeDiff = max(cell2time) - min(cell1time);
        timeDiff_true = time1(max(cell2time)) - time1(min(cell1time));
        dist_time(j) = timeDiff_true;                             % Distance in time
        maniTime = manifold3d([min(cell1time):max(cell2time)],:);
        dist_mani(j) = sum(sqrt(sum(diff(maniTime).^2,2)));  %Distance on manifold
        
        % for comparison, take a snippet of the same duration from a random
        % other trial on the manifold.
        for sh = 1:numShuf
            randtial = 0;
            while sum(trialn == randtial) <= (timeDiff + t0)
                randtial = datasample(unique(trialn),1);
            end
            if t0>0
                mi = min(find(trialn == randtial)) + randi(t0);
            else
                mi = min(find(trialn == randtial));
            end
            ma = mi + timeDiff;
            maniTime_rnd = manifold3d([mi:ma],:);
            dist_mani_rnd(j,sh) = sum(sqrt(sum(diff(maniTime_rnd).^2,2)));
        end
    end
    c = corrcoef(dist_mani, dist_time);
    corrs_per_doublet(i) = c(2,1);
    
    shuff_tmp = zeros(numShuf,1);
    for sh = 1:numShuf
        d = corrcoef(dist_mani_rnd(:,sh), dist_time);
        corrs_per_doublet_rnd(i,sh) = d(2,1);
    end
    
    dist_mani_all{i} = dist_mani;
    dist_time_all{i} = dist_time;
    dist_mani_rnd_all{i} = dist_mani_rnd;
end
close(f);
real_corrs = corrs_per_doublet;
shuffled_corrs = corrs_per_doublet_rnd;

outputDoubletManiCorrs.real_corrs = real_corrs;
outputDoubletManiCorrs.shuffled_corrs = shuffled_corrs;
outputDoubletManiCorrs.dist_mani_all = dist_mani_all;
outputDoubletManiCorrs.dist_time_all = dist_time_all;
outputDoubletManiCorrs.dist_mani_rnd_all = dist_mani_rnd_all;
outputDoubletManiCorrs.dist_mani_all_vertcat = vertcat(dist_mani_all{:});
outputDoubletManiCorrs.dist_time_all_vertcat = vertcat(dist_time_all{:});
outputDoubletManiCorrs.argins = argins;

