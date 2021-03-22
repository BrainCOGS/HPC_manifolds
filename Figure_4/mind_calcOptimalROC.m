function outputOptimalROC = mind_calcOptimalROC(fname, fname_mani, outputDoublets, numDim)
% numDim only applies to manifold, it's not used for other variables

argins.fname = fname;
argins.fname_mani = fname_mani;
argins.numDim = numDim;

load(fname_mani);

config.sd_range = 1:5:100;

nic_output = extractVariables('all', 2, 'keepTrials', 1, 0, 0, 0, [0 0], fname,'none','towers',1,1);
ROIactivities = nic_output.ROIactivities;
behavioralVariables = nic_output.behavioralVariables;

numSD = length(config.sd_range);
sd_range = config.sd_range;
totalROC_ALL = zeros(numSD,2);

for i=1:numSD
    curSD = sd_range(i);
    outputROC_ALL{i} = mind_calcDoubletROC(fname, outMind, outputDoublets, outputDoublets.saveAll_basics.sigAbove3, curSD,0, [], numDim);
    totalROC_ALL(i,:) = [outputROC_ALL{i}.totalTPR outputROC_ALL{i}.totalFPR];
    disp(['Completed ' num2str(i) ' of ' num2str(numSD) ' SDs.']);
end

totalOptimal = [ones(40,1) zeros(40,1)];
distancesROC = diag(pdist2(totalOptimal, totalROC_ALL,'euclidean'));
[~, m2] = min(distancesROC);



%% Save variables

outputOptimalROC.argins = argins;
outputOptimalROC.config = config;
outputOptimalROC.bestSD = config.sd_range(m2);
outputOptimalROC.ValuesROC = totalROC_ALL(m2,:);
outputOptimalROC.distancesROC = distancesROC;
outputOptimalROC.totalROC_ALL = totalROC_ALL;
outputOptimalROC.outputROC_ALL = outputROC_ALL;