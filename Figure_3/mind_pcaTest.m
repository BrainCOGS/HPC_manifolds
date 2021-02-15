function mean_trials_animals = mind_pcaTest(fnameStruct)

rng(1);

data_pca = NaN(7,10,150);
for animal = 1:length(fnameStruct)

    fname = fnameStruct(animal).fname;
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);

    ROIactivities = nic_output.ROIactivities;
    behavioralVariables = nic_output.behavioralVariables;
    trials = behavioralVariables.Trial;
    trialn = unique(trials);
    numtrials = length(trialn);
    CV = generateCrossValSet_v2(behavioralVariables, numtrials);
    CVtrials = randi(numtrials,10,1);
    
    for j=1:length(CVtrials)
        curTrial = CVtrials(j);
        trainingdata = ROIactivities(CV(curTrial).trainLocations,:);
        testdata     = ROIactivities(CV(curTrial).testLocations,:);
        [coeff,score,latent,tsquared,explained,mu1] = pca(trainingdata);
        
        for n_pc = 2:150
            score_testdata = testdata*coeff;
            rec_testdata = score_testdata(:,1:n_pc)*coeff(:,1:n_pc)';
            cc_test = corrcoef(rec_testdata(:), testdata(:));
            data_pca(animal,j,n_pc) = cc_test(2,1);
        end
    end
end


mean_trials = squeeze(mean(data_pca, 2));
mean_trials_animals = squeeze(mean(mean_trials, 1));

