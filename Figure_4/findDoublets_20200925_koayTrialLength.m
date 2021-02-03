function outputDoublets = findDoublets_20200925_koayTrialLength(stdcutoff, numDoublets, shuffletoggle, fname, preprocessparam, inputRNG)
% If detection type is 1, use edward version, if anything else, use sue-ann
% version
% If shuffltoggle is 1, shuffle the trialIDs
% if numDoublets is 0, it's going to give all doublets greater than 5
% instances
% preprocessparam is [3 0] or [5 2], or something like that

load(fname);
config.stdcutoff = stdcutoff;
argins.rng       = inputRNG;

rng(argins.rng);
%
% nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [5 2], fname,'none','towers',1,1);
nic_output = extractVariables('all', 2, 'keepTrials', 1, 0, 0, config.stdcutoff, preprocessparam, fname,'none','towers',1,1);

behavioralVariables = nic_output.behavioralVariables;
ROIactivities = nic_output.ROIactivities;
Evidence      = behavioralVariables.Evidence;

trialn = behavioralVariables.Trial;
[data_trial_nrs,~,ic] = unique(trialn);
ntrials= max(ic);

% dataDFF_trials = reshape(ROIactivities,[301,ntrials,size(ROIactivities,2)]);




% Number of times to shuffle the rois within doublets
shuffletimes = 100;
config.shufflestimes = shuffletimes;
config.stdcutoff     = stdcutoff;
config.numDoublets   = numDoublets;
config.shuffletoggle = shuffletoggle;
config.preprocessparam = preprocessparam;

% Find isGood and mainTrials to analyze
aTrials = data_trial_nrs;
score_subtrials = score.trial(aTrials);

% Use only the good trials
alltrials_choice = [score.trial(:).choice];
alltrials_choice = alltrials_choice(aTrials);
% alltrials_choice = [alltrials_choice(2:end) alltrials_choice(1)];
% alltrials_choice = alltrials_choice(randperm(length(alltrials_choice)));
lefttrials_choice = find(alltrials_choice==1);
righttrials_choice = find(alltrials_choice==2);
alltrials_trialType = [score.trial(:).trialType];
alltrials_trialType = alltrials_trialType(aTrials);
% alltrials_trialType = [alltrials_trialType(2:end) alltrials_trialType(1)];
% alltrials_trialType = alltrials_trialType(randperm(length(alltrials_trialType)));
lefttrials_trialType = find(alltrials_trialType==1);
righttrials_trialType = find(alltrials_trialType==2);

saveAll_basics.alltrials_choice = alltrials_choice;
saveAll_basics.lefttrials_choice = lefttrials_choice;
saveAll_basics.righttrials_choice = righttrials_choice;
saveAll_basics.alltrials_trialType = alltrials_trialType;
saveAll_basics.lefttrials_trialType = lefttrials_trialType;
saveAll_basics.righttrials_trialType = righttrials_trialType;
saveAll_basics.ROIactivities         = ROIactivities;

% Set the random number generator

saveAll_basics.stdcutoff = stdcutoff;


% for i=1:length(data_trial_nrs)
%     temp3a{i}  = ROIactivities(trialn==data_trial_nrs(i),:);
%     %temp3a_E{i}= Evidence(trialn==data_trial_nrs(i));
% %     temp4{1,i} = [];
% %     for j=1:size(temp3a{i},1)
% %         
% %         if ~isempty(find(temp3a{i}(j,:)==1))
% %             temp4{1,i} = [temp4{1,i} find(temp3a{i}(j,:)==1)];
% %         end
% %     end
% end

temp3a = computeXValidatedSequence_Edward(fname, stdcutoff, preprocessparam);


if shuffletoggle==1
    temp3a = seqanal_shuffleID(temp3a,lefttrials_choice, righttrials_choice);
end

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


%% Calculates the number of trials
bothinorderSaveLength = cellfun('length',bothinorderSave);

%% Find the trials where single cells are active

trialsROIActive = cell(length(score.roi),1);
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
saveAll_basics.aTrials = aTrials;
saveAll_basics.builtDoublets = builtDoublets;

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


% Nsamples = size(shufMean(:));
% alpha = 0.01;
% thresholdSigma = sqrt(2)*erfinv( (Nsamples - 2*alpha) / Nsamples);
% shufThreshold = shufMean+(thresholdSigma.*shufSTD);
shufThreshold = shufMean+(2.*shufSTD);

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
clear best1;
best1(1,:) = seq_row;
best1(2,:) = seq_col;
b(find(best1(1,:)==best1(2,:)))=[];
best1(:,find(best1(1,:)==best1(2,:)))=[];


bestall = best1;

% better5 is the number of doublets with more than 5 instances
% better5 = length(find(b>=5));
% best1=best1(:,1:better5);
% saveAll_basics.better5 = better5;
% better3 = length(find(b>=3));
% best1 = best1(:,1:better3);
% b     = b(1:better3);

% Get rid of nonsignificant doublets
best1 = best1(:,find(b));


% caps number of doublets if numDoublets is set
if numDoublets~=0
    if length(best1)>numDoublets
        best1=bestall(:,1:numDoublets);
    end
end

clear saveOrdered
clear saveLength


for i=1:size(best1,2)
    bothinorder1 = bothinorderSave{best1(1,i), best1(2,i)};
    bothinorder2 = alltrials_choice(bothinorder1);
    bothinorder2_trialType = alltrials_trialType(bothinorder1);
    saveOrdered(i) = length(find(bothinorder2==1))/length(bothinorder2);
    saveLength(i) = length(bothinorder2);
    saveOrdered_trialType(i) = length(find(bothinorder2_trialType==1))/length(bothinorder2_trialType);
    saveLength_trialType(i) = length(bothinorder2_trialType);
    saveTrials{i} = bothinorder1;
    saveTrials2{i} = bothinorder2;
end


% saveAll_trials = struct('activeROIs_hard', temp4a(1,:), 'activeROIs_soft', temp4(1,:));
saveAll_trials = struct('activeROIs_soft', temp4(1,:));
saveAll_doublets = struct('first_cell', num2cell(best1(1,:)), ...
    'second_cell', num2cell(best1(2,:)), ...
    'cells', num2cell(best1,1), ...
    'trials_appear', saveTrials, ...
    'trials_choices', saveTrials2, ...
    'number_doublets', num2cell(saveLength), ...
    'prediction',num2cell(saveOrdered));


% This loop takes about 4 minutes for E22

tic
wb=waitbar(0,'Calculating Significance of each Doublet Prediction');
for i=1:saveAll_doublets(1).number_doublets
    
    saveScramble = zeros(1000,1);
    for h=1:1000
        randlist = randsample(alltrials_choice, i);
        save1 = length(find(randlist==1))/i;
        saveScramble(h) = save1;
    end
    
    thres1(i) = std(saveScramble)*2+mean(saveScramble);
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



%% New shuffle where spike times are the same, but trial IDs are shuffled

[~, ~, BafterA_r_shf, BafterA_l_shf] = shuffleDoublets_v2( shuffletimes, alltrials_choice, lefttrials_choice, righttrials_choice, temp3a);

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


%% Calculate similarity of trials in which doublets show up versus random selection of the same # of L and R trials
% seqanal_trialdistance was changed, and doesn't work this with code now

% trialScores = seqanal_trialscore(score_subtrials, linspace(0,200,5));
%
% for i=1:length(saveAll_doublets)
%
%     exp_trials = saveAll_doublets(i).trials_appear;
%     exp_distance = seqanal_trialdistance(exp_trials, trialScores);
%
%     num_left = length(find(saveAll_doublets(i).trials_choices==1));
%     num_right = length(find(saveAll_doublets(i).trials_choices==2));
%
%     for j=1:30
%
%         rand_L = randperm(length(lefttrials_choice),num_left);
%         rand_R = randperm(length(righttrials_choice),num_right);
%
%         rand_L1 = lefttrials_choice(rand_L);
%         rand_R1 = righttrials_choice(rand_R);
%
%         rand_trials = [rand_L1 rand_R1];
%
%         rand_score(j) = seqanal_trialdistance(rand_trials, trialScores);
%
%     end
%
%     exp_score(i) = exp_distance;
%     rand_score_mean(i) = mean(rand_score);
%     rand_score_std(i) = std(rand_score);
%     score_sig(i) = exp_distance<mean(rand_score)-(2*std(rand_score));
%
% end
%
% exp_scoreCell = num2cell(exp_score);
% [saveAll_doublets.trials_distance] = exp_scoreCell{:};
% rand_score_meanCell = num2cell(rand_score_mean);
% [saveAll_doublets.random_distance_mean] = rand_score_meanCell{:};
% rand_score_stdCell = num2cell(rand_score_std);
% [saveAll_doublets.random_distance_std] = rand_score_stdCell{:};
% score_sigCell = num2cell(score_sig);
% [saveAll_doublets.score_sig] = score_sigCell{:};


%% Can use this code to get significant left and right doublets
saveAll_doublets_Sig = saveAll_doublets;
saveAll_doublets_Sig(find([saveAll_doublets(:).sig_predict]==0))=[];
saveAll_doublets_great5 = saveAll_doublets;
saveAll_doublets_great5(find([saveAll_doublets(:).number_doublets]<5))=[];

saveAll_doublets_NotSig = saveAll_doublets;
saveAll_doublets_NotSig(find([saveAll_doublets(:).sig_predict]==1))=[];
saveAll_doublets_NotSig(find([saveAll_doublets_NotSig(:).number_doublets]<5))=[];

saveAll_doublets_Sig_Left = saveAll_doublets_Sig;
saveAll_doublets_Sig_Left(find([saveAll_doublets_Sig_Left(:).direction]==1))=[];
saveAll_doublets_Sig_Right = saveAll_doublets_Sig;
saveAll_doublets_Sig_Right(find([saveAll_doublets_Sig_Right(:).direction]==0))=[];


saveAll_basics.number_doublets = [saveAll_doublets_great5.number_doublets];
saveAll_basics.shuffle_length = [saveAll_doublets_great5.shuffle_length];
saveAll_basics.number_sig = length(saveAll_doublets_Sig);

%% Make the outputfile

outputDoublets.saveAll_doublets = saveAll_doublets;
outputDoublets.saveAll_trials = saveAll_trials;
outputDoublets.saveAll_preprocess = saveAll_preprocess;
outputDoublets.saveAll_basics = saveAll_basics;
outputDoublets.bothinorderSave = bothinorderSave;
outputDoublets.score_subtrials = score_subtrials;
outputDoublets.saveAll_doublets_Sig = saveAll_doublets_Sig;
outputDoublets.saveAll_doublets_NotSig = saveAll_doublets_NotSig;
outputDoublets.saveAll_doublets_great5 = saveAll_doublets_great5;
outputDoublets.saveAll_doublets_Sig_Left = saveAll_doublets_Sig_Left;
outputDoublets.saveAll_doublets_Sig_Right = saveAll_doublets_Sig_Right;
outputDoublets.config = config;

findDoublets_20190509_plotter(outputDoublets)
