function outputDoublets = findDoublets(stdcutoff, fname, preprocessparam, inputRNG)
% If detection type is 1, use edward version, if anything else, use sue-ann
% version
% If shuffltoggle is 1, shuffle the trialIDs
% if numDoublets is 0, it's going to give all doublets greater than 5
% instances
% preprocessparam is [3 0] or [5 2], or something like that

%% Save inputs and set up config

argins.stdcutoff 	   = stdcutoff;
argins.fname           = fname;
argins.preprocessparam = preprocessparam;
argins.inputRNG        = inputRNG;

rng(inputRNG);

% Number of times to shuffle the rois within doublets
shuffletimes = 100;
config.shufflestimes = shuffletimes;

outputDoublets.argins = argins;
outputDoublets.config = config;

%% Set up behavioral variables

load(fname);
nic_output = extractVariables('all', 2, 'keepTrials', 1, 0, 0, stdcutoff, preprocessparam, fname,'none','towers',1,1);
behavioralVariables = nic_output.behavioralVariables;
ROIactivities = nic_output.ROIactivities;
trialn = behavioralVariables.Trial;
data_trial_nrs = unique(trialn);

% Save this for the findTriplets code
score_subtrials = score.trial(data_trial_nrs);
outputDoublets.score_subtrials = score_subtrials;

alltrials_choice = [score.trial(:).choice];
alltrials_choice = alltrials_choice(data_trial_nrs);
lefttrials_choice = find(alltrials_choice==1);
righttrials_choice = find(alltrials_choice==2);
alltrials_trialType = [score.trial(:).trialType];
alltrials_trialType = alltrials_trialType(data_trial_nrs);

saveAll_basics.alltrials_choice = alltrials_choice;
saveAll_basics.lefttrials_choice = lefttrials_choice;
saveAll_basics.righttrials_choice = righttrials_choice;
saveAll_basics.alltrials_trialType = alltrials_trialType;
saveAll_basics.ROIactivities         = ROIactivities;
saveAll_basics.data_trialnrs = data_trial_nrs;

clear score;


%% Find doublets
for i=1:length(data_trial_nrs)
    temp3a{i}  = ROIactivities(trialn==data_trial_nrs(i),:);
end

% Be careful that these is multiple ROIs are active at same time bin, it
% just put them in descending order
for i=1:length(data_trial_nrs)  
    temp4{1,i} = [];
    for j=1:size(temp3a{i},1)
        
        if ~isempty(find(temp3a{i}(j,:)==1))
            temp4{1,i} = [temp4{1,i} find(temp3a{i}(j,:)==1)];
        end
    end
end

saveAll_preprocess = struct('digitized_eventStart', temp3a, ...
                            'orderedEvents',temp4(1,:));


%% Find the doublets, computes trials in which the pair of ROIs shows up

[bothinorderSave, builtDoublets] = buildDoublets(temp3a);
saveAll_basics.builtDoublets = builtDoublets;

%% Calculates the number of trials

bothinorderSaveLength = cellfun('length',bothinorderSave);


%% Find the trials where single cells are active (used in findTriplets)

trialsROIActive = cell(size(ROIactivities,2),1);
for i=1:length(alltrials_choice)
    
    curTrial = saveAll_preprocess(i).orderedEvents;
    
    for j=1:length(curTrial)
        if isempty(trialsROIActive{curTrial(j)})
            trialsROIActive{curTrial(j)} = i;
        elseif ~ismember(i, trialsROIActive{curTrial(j)})
            trialsROIActive{curTrial(j)} = [trialsROIActive{curTrial(j)} i];
        end
    end
    
end
saveAll_basics.trialsROIActive = trialsROIActive;


%% Shuffle the activity for each ROI independently in each trial

bothinorderSaveLength_shuf = zeros([100 size(bothinorderSave)]);
for k=1:100
    for i=1:length(temp3a)
        
        trialData = temp3a{i};
        
        shufData = zeros(size(trialData));
        for j=1:size(trialData,2)
            
            roiData = trialData(:,j);
            shiftby = 1+randi(length(roiData)-2);
            shufData(:,j) = circshift(roiData,shiftby);
        end
        
        temp3a_shuf{i} = shufData;
        
    end
    
    [bothinorderSave_shuf, ~] = buildDoublets(temp3a_shuf);
    bothinorderSaveLength_shuf(k,:,:) = cellfun('length',bothinorderSave_shuf);
end

shufMean = squeeze(mean(bothinorderSaveLength_shuf,1));
shufSTD  = squeeze(std(bothinorderSaveLength_shuf,1));
shufThreshold = shufMean+(2.*shufSTD);

% Get rid of all doublets that were not greater than the shuffled threshold
% and also did not appear more than 3 times
sigDoublets = bothinorderSaveLength>shufThreshold;
bothinorderSaveLength_3andAbove = bothinorderSaveLength>=3;
sigAbove3 = bothinorderSaveLength_3andAbove & sigDoublets;

realLength = bothinorderSaveLength;
realLength(~sigAbove3)=0;
bothinorderSaveLength = realLength;

saveAll_basics.sigDoublets = sigDoublets;
saveAll_basics.sigAbove3   = sigAbove3;

%% Find the top doublets

[b, ix] = sort(bothinorderSaveLength(:),'descend');
[seq_row, seq_col] = ind2sub(size(bothinorderSaveLength),ix);
best1(1,:) = seq_row;
best1(2,:) = seq_col;

% Get rid of doublets where it's the same cell as first and second cell
sameLoc = find(best1(1,:)==best1(2,:));
b(sameLoc)=[];
best1(:,sameLoc)=[];

% Get rid of nonsignificant and doublets appearing fewer then 3 times
best1 = best1(:,find(b));

%% Get important doublet statistics

for i=1:size(best1,2)
    bothinorder1 = bothinorderSave{best1(1,i), best1(2,i)};
    bothinorder2 = alltrials_choice(bothinorder1);
    bothinorder2_trialType = alltrials_trialType(bothinorder1);
    saveOrdered(i) = length(find(bothinorder2==1))/length(bothinorder2);
    saveLength(i) = length(bothinorder2);
    saveTrials{i} = bothinorder1;
    saveTrials2{i} = bothinorder2;
end

saveAll_trials = struct('activeROIs_soft', temp4);
saveAll_doublets = struct('first_cell', num2cell(best1(1,:)), ...
    'second_cell', num2cell(best1(2,:)), ...
    'cells', num2cell(best1,1), ...
    'trials_appear', saveTrials, ...
    'trials_choices', saveTrials2, ...
    'number_doublets', num2cell(saveLength), ...
    'prediction',num2cell(saveOrdered));


%% Calculate the significance of each doublet's prediction

tic
wb=waitbar(0,'Calculating Significance of each Doublet Prediction');

% Find the thresholds for predictiveness of a doublet that occurs i times
for i=1:saveAll_doublets(1).number_doublets
    
    saveScramble = zeros(1000,1);
    for h=1:1000
        randlist = randsample(alltrials_choice, i);
        save1 = length(find(randlist==1))/i;
        saveScramble(h) = save1;
    end
    
    thres1(i) = mean(saveScramble) + (std(saveScramble)*2);
    thres2(i) = mean(saveScramble) - (std(saveScramble)*2);
end

for j=1:size(best1,2)
    
    if saveAll_doublets(j).prediction>thres1(saveAll_doublets(j).number_doublets) || saveAll_doublets(j).prediction<thres2(saveAll_doublets(j).number_doublets)
        sig1(j)=1;
    else
        sig1(j)=0;
    end
    
    currentFirst = saveAll_doublets(j).cells(1);
    currentSecond = saveAll_doublets(j).cells(2);
    onlyFirst=[];
    onlySecond=[];
    allFirst=[];
    allSecond=[];
    
    for i=1:length(saveAll_trials)
        
        if ~isempty(find(saveAll_trials(i).activeROIs_soft==currentFirst)) & isempty(find(saveAll_trials(i).activeROIs_soft==currentSecond))
            onlyFirst = [onlyFirst i];
        end
        if isempty(find(saveAll_trials(i).activeROIs_soft==currentFirst)) & ~isempty(find(saveAll_trials(i).activeROIs_soft==currentSecond))
            onlySecond = [onlySecond i];
        end
        
        if ~isempty(find(saveAll_trials(i).activeROIs_soft==currentFirst))
            allFirst = [allFirst i];
        end
        if isempty(find(saveAll_trials(i).activeROIs_soft==currentFirst))
            allSecond = [allSecond i];
        end
        
    end
    
    onlyFirst1{j} = onlyFirst;
    onlySecond1{j} = onlySecond;
    
    onlyFirst2{j} = alltrials_choice(onlyFirst);
    onlySecond2{j} = alltrials_choice(onlySecond);
    
    onlyFirst3{j} = length(find(onlyFirst2{j}==1))/length(onlyFirst2{j});
    onlySecond3{j} = length(find(onlySecond2{j}==1))/length(onlySecond2{j});
    
    
    allFirst1{j} = allFirst;
    allSecond1{j} = allSecond;
    
    allFirst2{j} = alltrials_choice(allFirst);
    allSecond2{j} = alltrials_choice(allSecond);
    
    allFirst3{j} = length(find(allFirst2{j}==1))/length(allFirst2{j});
    allSecond3{j} = length(find(allSecond2{j}==1))/length(allSecond2{j});
    
    waitbar(j/size(best1,2));
end
close(wb);
disp('Calculating significance');
toc

[saveAll_doublets.only_first] = onlyFirst1{:};
[saveAll_doublets.only_second] = onlySecond1{:};
[saveAll_doublets.only_first_choice] = onlyFirst2{:};
[saveAll_doublets.only_second_choice] = onlySecond2{:};
[saveAll_doublets.only_first_predict] = onlyFirst3{:};
[saveAll_doublets.only_second_predict] = onlySecond3{:};

[saveAll_doublets.all_first] = allFirst1{:};
[saveAll_doublets.all_second] = allSecond1{:};
[saveAll_doublets.all_first_choice] = allFirst2{:};
[saveAll_doublets.all_second_choice] = allSecond2{:};
[saveAll_doublets.all_first_predict] = allFirst3{:};
[saveAll_doublets.all_second_predict] = allSecond3{:};


%% Shuffle where spike times are the same, but trial IDs are shuffled

[~, ~, BafterA_r_shf, BafterA_l_shf] = shuffleDoublets_v2(shuffletimes, alltrials_choice, lefttrials_choice, righttrials_choice, temp3a);

for i=1:length(best1)
    
    lefts(i) = sum([BafterA_l_shf(best1(1,i),best1(2,i),:)]);
    rights(i) = sum([BafterA_r_shf(best1(1,i),best1(2,i),:)]);
    
    together1 = [BafterA_l_shf(best1(1,i),best1(2,i),:)]+[BafterA_r_shf(best1(1,i),best1(2,i),:)];
    together1_reverse = [BafterA_l_shf(best1(2,i),best1(1,i),:)]+[BafterA_r_shf(best1(2,i),best1(1,i),:)];
    means1_diff(i) = mean(together1 - together1_reverse);
    means1_reverse(i) = mean(together1_reverse);
    means1(i) = mean(together1);
    stds1(i) = std(together1);
    
    tots(i) = lefts(i)+rights(i);
    fracs(i) = lefts(i)/tots(i);
    
    if saveAll_doublets(i).trials_appear>(means1(i) + (2*stds1(i)))
        sigNum1(j)=1;
    else
        sigNum1(j)=0;
    end
    
end

savedShuffles3cell = num2cell(fracs);
[saveAll_doublets.doublet_shuffle] = savedShuffles3cell{:};
sig1cell = num2cell(sig1);
[saveAll_doublets.sig_predict] = sig1cell{:};
sigNum1cell = num2cell(sigNum1);
[saveAll_doublets.sig_num] = sigNum1cell{:};
shuffleLength1 = num2cell(tots./shuffletimes);
[saveAll_doublets.shuffle_length] = shuffleLength1{:};

mean_shuf = num2cell(means1);
[saveAll_doublets.mean_shuffle_length] = mean_shuf{:};
mean_reverse_shuf = num2cell(means1_reverse);
[saveAll_doublets.mean_shuffle_reverse_length] = mean_reverse_shuf{:};
std_shuf = num2cell(stds1);
[saveAll_doublets.std_shuffle_length] = std_shuf{:};
asym_diff_shuf = num2cell(means1_diff);
[saveAll_doublets.mean_asymmetry_diff_length_shuf] = asym_diff_shuf{:};

%% Set the direction of each doublet, 1 means right-predicting

for i=1:length(saveAll_doublets)
    
    if saveAll_doublets(i).prediction<.5
        saveAll_doublets(i).direction=1;
    else
        saveAll_doublets(i).direction=0;
    end
    
end


%% Asymmetric-ness of the doublets
for i=1:length(saveAll_doublets)
    
    num_asym = cell2mat(bothinorderSave(saveAll_doublets(i).cells(2), saveAll_doublets(i).cells(1)));
    num_asym_length = length(num_asym);
    
    asym{i} = num_asym;
    asym_l(i) = num_asym_length;
    
end

asym_l_cell = num2cell(asym_l);
[saveAll_doublets.asymmetry_length] = asym_l_cell{:};
[saveAll_doublets.asymmetry_trials] = asym{:};

asym_diff = num2cell(saveLength-asym_l);
[saveAll_doublets.asymmetry_length_diff] = asym_diff{:};


%% Can use this code to get significant left and right doublets
saveAll_doublets_Sig = saveAll_doublets;
saveAll_doublets_Sig(find([saveAll_doublets(:).sig_predict]==0))=[];

saveAll_doublets_NotSig = saveAll_doublets;
saveAll_doublets_NotSig(find([saveAll_doublets(:).sig_predict]==1))=[];

saveAll_doublets_Sig_Left = saveAll_doublets_Sig;
saveAll_doublets_Sig_Left(find([saveAll_doublets_Sig_Left(:).direction]==1))=[];
saveAll_doublets_Sig_Right = saveAll_doublets_Sig;
saveAll_doublets_Sig_Right(find([saveAll_doublets_Sig_Right(:).direction]==0))=[];

saveAll_basics.number_sig = length(saveAll_doublets_Sig);

%% Make the outputfile

outputDoublets.saveAll_doublets = saveAll_doublets;
outputDoublets.saveAll_trials = saveAll_trials;
outputDoublets.saveAll_preprocess = saveAll_preprocess;
outputDoublets.saveAll_basics = saveAll_basics;
outputDoublets.bothinorderSave = bothinorderSave;
outputDoublets.saveAll_doublets_Sig = saveAll_doublets_Sig;
outputDoublets.saveAll_doublets_NotSig = saveAll_doublets_NotSig;
outputDoublets.saveAll_doublets_Sig_Left = saveAll_doublets_Sig_Left;
outputDoublets.saveAll_doublets_Sig_Right = saveAll_doublets_Sig_Right;

