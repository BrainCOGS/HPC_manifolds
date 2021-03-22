function outputROC = mind_calcDoubletROC(fname, outMind, outputDoublets, sigAbove3, stdCutoff, toggleShuffle, full_reconstruction, numDim)
% if toggleShuffle ==1, shuffle the data within trial
% numDim is the dimensions of the embedded manifold
argins.fname = fname;
argins.stdCutoff = stdCutoff;
argins.toggleShuffle = toggleShuffle;
argins.full_reconstruction = full_reconstruction;
argins.numDim = numDim;

nic_output = extractVariables('all', 2, 'keepTrials', 1, 0, 0, 0, [0 0], fname,'none','towers',1,1);
ROIactivities = nic_output.ROIactivities;
behavioralVariables = nic_output.behavioralVariables;
trialn = behavioralVariables.Trial;
[data_trial_nrs,~,ic] = unique(trialn);
numTrials = length(data_trial_nrs);

if isempty(full_reconstruction)
    
    full_reconstruction = zeros(size(ROIactivities));
    
    sample = outMind.dataDFF;
    pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
    dat = outMind.dat;
    a = dat.allembed([outMind.dat.allembed.d]==numDim).m2f.map.transform( dat.allembed([outMind.dat.allembed.d]==numDim).f2m.map.transform(pca_coords) );
    b = dat.forestdat.pca.model.inverse_transform(a);
    
    full_reconstruction(outMind.Datarange,outMind.Neurons) = b;
end

noiseROI = zeros(1,size(full_reconstruction,2));
for i=1:size(full_reconstruction,2)
    noiseROI(i) = robustSTD(full_reconstruction(:,i));
end

bothinorderSave = outputDoublets.bothinorderSave;
bothinorderSaveLength = cellfun('length',bothinorderSave);
temp3a             = {outputDoublets.saveAll_preprocess.digitized_eventStart};

lefttrials_choice  = outputDoublets.saveAll_basics.lefttrials_choice;
righttrials_choice = outputDoublets.saveAll_basics.righttrials_choice;

%% Calculate the statistics

DFFthresholded = full_reconstruction > (argins.stdCutoff * noiseROI);

% take diff to find event starts where diff is 1
DFFdiff = diff(DFFthresholded);
DFFdiff(DFFdiff == -1) = 0;
% diff loses a row, so use the first row of DFFthresholded as the
% first row.
DFF = [DFFthresholded(1,:); DFFdiff];

for i=1:length(data_trial_nrs)
    temp3a_mani{i}  = DFF(trialn==data_trial_nrs(i),:);
end

if toggleShuffle==1
    temp3a_mani = seqanal_shuffleID(temp3a_mani,lefttrials_choice, righttrials_choice);
end

[bothinorderSave_mani, ~] = buildDoublets(temp3a_mani);

bothinorderSaveLength_mani = cellfun('length',bothinorderSave_mani);

truepositive_falsenegative = cellfun(@ismember,bothinorderSave, bothinorderSave_mani,'UniformOutput',false);
falsepositiveIs0           = cellfun(@ismember,bothinorderSave_mani, bothinorderSave,'UniformOutput',false);

numTP = cellfun(@sum,truepositive_falsenegative);
numFN = bothinorderSaveLength - numTP;
numTP2 = cellfun(@sum,falsepositiveIs0);
numFP = bothinorderSaveLength_mani - numTP2;

totalTrials = ones(size(bothinorderSave));
totalTrials = totalTrials.*numTrials;
numTN = totalTrials - (numTP+numFN+numFP);

FPR = numFP./(numFP+numTN);
TPR = numTP./(numTP+numFN);

ACC = (numTP+numTN)./(numTP+numTN+numFP+numFN);
PRE = numTP./(numTP+numFP);

%% Save variables

FPR_doublet = FPR(sigAbove3);
TPR_doublet = TPR(sigAbove3);
numTP_doublet = numTP(sigAbove3);
numFP_doublet = numFP(sigAbove3);
numFN_doublet = numFN(sigAbove3);
numTN_doublet = numTN(sigAbove3);
ACC_doublet   = ACC(sigAbove3);
PRE_doublet   = PRE(sigAbove3);

meanTPR = mean(TPR_doublet);
meanFPR = mean(FPR_doublet);
meanTP  = mean(numTP_doublet);
meanFP  = mean(numFP_doublet);
meanTN  = mean(numTN_doublet);
meanFN  = mean(numFN_doublet);
meanACC = mean(ACC_doublet);
meanPRE = nanmean(PRE_doublet);

totalTPR = sum(numTP_doublet)/(sum(numTP_doublet)+sum(numFN_doublet));
totalFPR = sum(numFP_doublet)/(sum(numFP_doublet)+sum(numTN_doublet));
totalPerDoublets = sum(numFN_doublet+numTP_doublet)/(length(numFN_doublet)*length(temp3a));

outputROC.argins = argins;
outputROC.sigAbove3   = sigAbove3;
outputROC.numTP       = numTP;
outputROC.numFP       = numFP;
outputROC.numTN       = numTN;
outputROC.numFN       = numFN;
outputROC.FPR         = FPR;
outputROC.TPR         = TPR;
outputROC.ACC         = ACC;
outputROC.PRE         = PRE;
outputROC.meanTPR = meanTPR;
outputROC.meanFPR = meanFPR;
outputROC.meanTP  = meanTP;
outputROC.meanFP  = meanFP;
outputROC.meanTN  = meanTN;
outputROC.meanFN  = meanFN;
outputROC.meanACC = meanACC;
outputROC.meanPRE = meanPRE;
outputROC.TPR_doublet = TPR_doublet;
outputROC.FPR_doublet = FPR_doublet;
outputROC.numTP_doublet = numTP_doublet;
outputROC.numFP_doublet = numFP_doublet;
outputROC.numTN_doublet = numTN_doublet;
outputROC.numFN_doublet = numFN_doublet;
outputROC.ACC_doublet   = ACC_doublet;
outputROC.PRE_doublet   = PRE_doublet;
outputROC.totalTPR = totalTPR;
outputROC.totalFPR = totalFPR;
outputROC.totalPerDoublets = totalPerDoublets;
outputROC.DFF    = DFF;
outputROC.full_reconstruction = full_reconstruction;

