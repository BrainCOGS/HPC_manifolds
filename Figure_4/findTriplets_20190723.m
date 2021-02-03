function outputTriplets = findTriplets_20190723(outputDoublets)
rng(1);
saveAll_preprocess = outputDoublets.saveAll_preprocess;
saveAll_doublets   = outputDoublets.saveAll_doublets;
saveAll_basics     = outputDoublets.saveAll_basics;
score_subtrials    = outputDoublets.score_subtrials;
temp3a             = {outputDoublets.saveAll_preprocess.digitized_eventStart};


%% Find the triplets

% [CafterBafterA, ball] = buildTriplets(temp3a,1);
% bothinorderSaveLengthT = cellfun('length',CafterBafterA);


% bothinorderSaveLengthT_shuf = false([100 size(CafterBafterA)]);
% bothinorderSaveLengthT_shuf = false([100 size(saveAll_basics.ROIactivities,2) size(saveAll_basics.ROIactivities,2) size(saveAll_basics.ROIactivities,2)]);
% for k=1:100
%     for i=1:length(temp3a)
%         
%         trialData = temp3a{i};
%         
%         shufData = zeros(size(trialData));
%         for j=1:size(trialData,2)
%             
%             roiData = trialData(:,j);
%             shiftby = 1+randi(length(roiData)-2);
%             shufData(:,j) = circshift(roiData,shiftby);
%         end
%         
%         temp3a_shuf{i} = shufData;
%         
%     end
%     
%     [~, ballT] = buildTriplets(temp3a_shuf,0);
%     bothinorderSaveLengthT_shuf(k,:,:,:) = sum(ballT,4);
%     clear ballT;
%     k
% end
% 
% shufMean = squeeze(mean(bothinorderSaveLengthT_shuf,1));
% shufSTD  = squeeze(std(bothinorderSaveLengthT_shuf,1));
% 
% % Nsamples = size(shufMean(:));
% % alpha = 0.01;
% % thresholdSigma = sqrt(2)*erfinv( (Nsamples - 2*alpha) / Nsamples);
% % shufThreshold = shufMean+(thresholdSigma.*shufSTD);
% 
% shufThreshold = shufMean+(2.*shufSTD);
% sigTriplets = bothinorderSaveLengthT>shufThreshold;
% 
% realLength = bothinorderSaveLengthT;
% realLength(~sigTriplets)=0;
% bothinorderSaveLengthT = realLength;
% [b, ix] = sort(bothinorderSaveLengthT(:),'descend');
% [seq1, seq2, seq3] = ind2sub(size(bothinorderSaveLengthT),ix);
% clear best1;
% best1(1,:) = seq1;
% best1(2,:) = seq2;
% best1(3,:) = seq3;
% b(find(best1(1,:)==best1(2,:) | best1(2,:)==best1(3,:) | best1(1,:)==best1(3,:)))=[];
% best1(:,find(best1(1,:)==best1(2,:) | best1(2,:)==best1(3,:) | best1(1,:)==best1(3,:)))=[];
% best1 = best1(:,find(b));
% sigb = b(find(b));






trialROIactive = saveAll_basics.trialsROIActive';
aTrials = saveAll_basics.aTrials;
alltrials_choice = saveAll_basics.alltrials_choice;
lefttrials_choice = saveAll_basics.lefttrials_choice;
righttrials_choice = saveAll_basics.righttrials_choice;
builtDoublets = saveAll_basics.builtDoublets;

% Stores whether each doublet appeared in each trial
trialsboolean = false(length(saveAll_doublets),length(score_subtrials));
for i=1:length(saveAll_doublets)
    trialsboolean(i,saveAll_doublets(i).trials_appear)=true;
end

% Stores whether each ROI appeared in each trial
roiboolean = false(length(saveAll_basics.trialsROIActive),length(score_subtrials));
for i=1:length(saveAll_basics.trialsROIActive)
    roiboolean(i,saveAll_basics.trialsROIActive{i})=true;
end


tic
wb=waitbar(0,'Looking for Triplets');
tCount = 1;
% tripletSave = cell(100000000,5);
for i=1:length(saveAll_doublets)
    cell1 = saveAll_doublets(i).first_cell;
    cell2 = saveAll_doublets(i).second_cell;
    posT1 = find([saveAll_doublets(:).first_cell]==cell2);
    
    for j=1:length(posT1)
        cell3 = saveAll_doublets(posT1(j)).second_cell;
        overlapTrials = trialsboolean(i,:).*trialsboolean(posT1(j),:);
        
        if sum(overlapTrials)>=3
            
            
            onlycell1 = (roiboolean(cell1,:)) & ~((roiboolean(cell2,:)+roiboolean(cell3,:)).*roiboolean(cell1,:));
            onlycell2 = (roiboolean(cell2,:)) & ~((roiboolean(cell1,:)+roiboolean(cell3,:)).*roiboolean(cell2,:));
            onlycell3 = (roiboolean(cell3,:)) & ~((roiboolean(cell1,:)+roiboolean(cell2,:)).*roiboolean(cell3,:));
            onlycell1and2 = squeeze(builtDoublets(cell1,cell2,:))' & ~(squeeze(builtDoublets(cell1,cell2,:))'.*roiboolean(cell3,:));
            onlycell2and3 = squeeze(builtDoublets(cell2,cell3,:))' & ~(squeeze(builtDoublets(cell2,cell3,:))'.*roiboolean(cell1,:));
            onlycell1and3 = squeeze(builtDoublets(cell1,cell3,:))' & ~(squeeze(builtDoublets(cell1,cell3,:))'.*roiboolean(cell2,:));
            
            %         onlycell1 = setdiff(trialROIactive{cell1},unique([trialROIactive{cell2} trialROIactive{cell3}]));
            %         onlycell2 = setdiff(trialROIactive{cell2},unique([trialROIactive{cell1} trialROIactive{cell3}]));
            %         onlycell3 = setdiff(trialROIactive{cell3},unique([trialROIactive{cell2} trialROIactive{cell1}]));
            %         onlycell1and2 = setdiff(bothinorderSave{cell1, cell2}, trialROIactive{cell3});
            %         onlycell2and3 = setdiff(bothinorderSave{cell2, cell3}, trialROIactive{cell1});
            %         onlycell1and3 = setdiff(bothinorderSave{cell1, cell3}, trialROIactive{cell2});
            
            tripletSave{tCount,1} = uint16(cell1);
            tripletSave{tCount,2} = uint16(cell2);
            tripletSave{tCount,3} = uint16(cell3);
            tripletSave{tCount,4} = uint16(sum(overlapTrials));
            tripletSave{tCount,5} = find(overlapTrials==1);
            tripletSave{tCount,6} = uint16(find(onlycell1));
            tripletSave{tCount,7} = uint16(find(onlycell2));
            tripletSave{tCount,8} = uint16(find(onlycell3));
            tripletSave{tCount,9} = uint16(find(onlycell1and2));
            tripletSave{tCount,10} = uint16(find(onlycell2and3));
            tripletSave{tCount,11} = uint16(find(onlycell1and3));
            tCount = tCount+1;
        end
    end
    waitbar(i/length(saveAll_doublets), wb, sprintf('%i', i));
end
close(wb);
disp('Finding Triplets');
toc

replist=[];
for i=1:size(tripletSave,1)
    rep1 = tripletSave{i,1};
    rep2 = tripletSave{i,3};
    if rep1==rep2
        replist = [replist i];
    end
    
end
tripletSave(replist,:)=[];
tripletSave = sortrows(tripletSave,-4);

%% Calculate the shuffles only for triplets

shufTimes = 100;
tripletMean = zeros(size(tripletSave,1),1);
tripletSTD = zeros(size(tripletSave,1),1);
for i=1:size(tripletSave,1)
    
    cell1 = tripletSave{i,1};
    cell2 = tripletSave{i,2};
    cell3 = tripletSave{i,3};
    listTrials = tripletSave{i,5};
    trialShufData = false(shufTimes,length(listTrials));
    for j=1:length(listTrials)
        
        curTrial = listTrials(j);
        trialActivity = temp3a{curTrial};
        roiActivity = trialActivity(:,[cell1 cell2 cell3]);
        
        for k=1:shufTimes
            
            shufData = zeros(length(roiActivity),3);
            for m=1:size(roiActivity,2)
                
                roiData = roiActivity(:,m);
                shiftby = 1+randi(length(roiData)-2);
                shufData(:,m) = circshift(roiData,shiftby);
            end
            cell1times = find(shufData(:,1));
            cell2times = find(shufData(:,2));
            cell3times = find(shufData(:,3));
            
            min1 = min(cell1times);
            max3 = max(cell3times);
            clear tempComp;
            clear tempComp2;
            tempComp = find(cell2times>min1);
            if sum(tempComp)>0
                tempComp2 = max3>cell2times(tempComp);
                
                if sum(tempComp2)>0
                    trialShufData(k,j)=true;
                end
            end
            
            
            
        end
        
        
    end
    
    tripletShufData = sum(trialShufData,2);
    tripletMean(i) = mean(tripletShufData);
    tripletSTD(i)  = std(tripletShufData);
    
end

shufThreshold = tripletMean+(2.*tripletSTD);

%% Take only triplets that pass the shuffled threshold

sigTrip = false(size(tripletSave,1),1);
for i=1:size(tripletSave,1) 
    curThres = shufThreshold(i);
    sigTrip(i) = tripletSave{i,4}>curThres;
end

tripletSave = tripletSave(sigTrip,:);

%% Take the top triplets

% better4 = length(find([tripletSave{:,4}]>=4));
% tripletSave = tripletSave(1:better4,:);

% using numTriplets
% tripletSave = tripletSave(1:numTriplets,:);

tic
clear saveOrdered_T
clear saveLength_T
clear saveChoices_T
for i=1:size(tripletSave,1)
    
    bothinorder1_T = cell2mat(tripletSave(i,5));
    bothinorder2_T = [alltrials_choice(bothinorder1_T)];
    saveOrdered_T(i) = length(find(bothinorder2_T==1))/length(bothinorder2_T);
    saveLength_T(i) = length(bothinorder2_T);
    saveChoices_T{i} = bothinorder2_T;
    
    onlycell1Choice = [alltrials_choice(cell2mat(tripletSave(i,6)))];
    onlycell2Choice = [alltrials_choice(cell2mat(tripletSave(i,7)))];
    onlycell3Choice = [alltrials_choice(cell2mat(tripletSave(i,8)))];
    
    %     onlycell1Trials = tripletSave(i,6);
    %     onlycell2Trials = tripletSave(i,7);
    %     onlycell3Trials = tripletSave(i,8);
    
    saveOrdered_cell1(i) = length(find(onlycell1Choice==1))/length(onlycell1Choice);
    saveOrdered_cell2(i) = length(find(onlycell2Choice==1))/length(onlycell2Choice);
    saveOrdered_cell3(i) = length(find(onlycell3Choice==1))/length(onlycell3Choice);
    
    saveLength_cell1(i) = length(onlycell1Choice);
    saveLength_cell2(i) = length(onlycell2Choice);
    saveLength_cell3(i) = length(onlycell3Choice);
    
    saveChoices_cell1{i} = onlycell1Choice;
    saveChoices_cell2{i} = onlycell2Choice;
    saveChoices_cell3{i} = onlycell3Choice;
    
    onlycell1and2Choice = [alltrials_choice(cell2mat(tripletSave(i,9)))];
    onlycell2and3Choice = [alltrials_choice(cell2mat(tripletSave(i,10)))];
    onlycell1and3Choice = [alltrials_choice(cell2mat(tripletSave(i,11)))];
    
    %     onlycell1and2Trials = tripletSave(i,9);
    %     onlycell2and3Trials = tripletSave(i,10);
    %     onlycell1and3Trials = tripletSave(i,11);
    
    saveOrdered_cell1and2(i) = length(find(onlycell1and2Choice==1))/length(onlycell1and2Choice);
    saveOrdered_cell2and3(i) = length(find(onlycell2and3Choice==1))/length(onlycell2and3Choice);
    saveOrdered_cell1and3(i) = length(find(onlycell1and3Choice==1))/length(onlycell1and3Choice);
    
    saveLength_cell1and2(i) = length(onlycell1and2Choice);
    saveLength_cell2and3(i) = length(onlycell2and3Choice);
    saveLength_cell1and3(i) = length(onlycell1and3Choice);
    
    saveChoices_cell1and2{i} = onlycell1and2Choice;
    saveChoices_cell2and3{i} = onlycell2and3Choice;
    saveChoices_cell1and3{i} = onlycell1and3Choice;
    
    settrials = [cell2mat(tripletSave(i,5)) cell2mat(tripletSave(i,9))];
    for j=1:100
        
        y = randsample(length(settrials),length(cell2mat(tripletSave(i,9))));
        
        subpredict(j) = sum(alltrials_choice(settrials(y))==1)/length(cell2mat(tripletSave(i,9)));
    end
    submean(i) = mean(subpredict);
    substd(i) = std(subpredict);
    
    if saveOrdered_T(i)>.5 & saveOrdered_cell1and2(i)<(submean(i)-(2*substd(i)))
        sub_sig(i) =  1;
    elseif saveOrdered_T(i)<.5 & saveOrdered_cell1and2(i)>(submean(i)+(2*substd(i)))
        sub_sig(i) =  1;
    else
        sub_sig(i) =  0;
    end
    
    
    
end
toc

%% Find significant over shuffle
tic
wb=waitbar(0,'Calculating Significance of each Triplet Prediction');
for i=1:double(tripletSave{1,4})
    
    saveScramble = zeros(1000,1);
    for h=1:1000
        randlist = randsample(alltrials_choice, i);
        save1 = length(find(randlist==1))/i;
        saveScramble(h) = save1;
    end
    
    thres1(i) = std(saveScramble)*2+mean(saveScramble);
    thres2(i) = mean(saveScramble) - (std(saveScramble)*2);
end

for j=1:size(tripletSave,1)
    
    
    if saveOrdered_T(j)>thres1(tripletSave{j,4}) || saveOrdered_T(j)<thres2(tripletSave{j,4})
        sig1T(1,j)=1;
    else
        sig1T(1,j)=0;
    end
    sig1T(2,j) = mean(saveScramble);
    sig1T(3,j) = std(saveScramble);
end
close(wb);
toc

% Create Structure for Triplets
saveAll_triplets = struct('first_cell', tripletSave(:,1), ...
    'second_cell', tripletSave(:,2),...
    'third_cell', tripletSave(:,3),...
    'cells', num2cell(cell2mat(tripletSave(:,1:3)),2), ...
    'trials_appear',tripletSave(:,5), ...
    'trials_choices',saveChoices_T',...
    'number_triplets', tripletSave(:,4), ...
    'prediction',num2cell(saveOrdered_T)',...
    'cell1only_predict',num2cell(saveOrdered_cell1)',...
    'cell2only_predict',num2cell(saveOrdered_cell2)',...
    'cell3only_predict',num2cell(saveOrdered_cell3)',...
    'cell1only_trials', tripletSave(:,6),...
    'cell2only_trials', tripletSave(:,7),...
    'cell3only_trials', tripletSave(:,8),...
    'cell1and2only_predict',num2cell(saveOrdered_cell1and2)',...
    'cell2and3only_predict',num2cell(saveOrdered_cell2and3)',...
    'cell1and3only_predict',num2cell(saveOrdered_cell1and3)',...
    'cell1and2only_trials', tripletSave(:,9),...
    'cell2and3only_trials', tripletSave(:,10),...
    'cell1and3only_trials', tripletSave(:,11),...
    'cell1only_length',num2cell(saveLength_cell1)',...
    'cell2only_length',num2cell(saveLength_cell2)',...
    'cell3only_length',num2cell(saveLength_cell3)',...
    'cell1and2only_length',num2cell(saveLength_cell1and2)',...
    'cell2and3only_length',num2cell(saveLength_cell2and3)',...
    'cell1and3only_length',num2cell(saveLength_cell1and3)',...
    'cell1and2only_shufflemean',num2cell(submean)',...
    'cell1and2only_shufflestd',num2cell(substd)',...
    'cell1and2only_shuffleSig',num2cell(sub_sig)');

%% New shuffle FOR TRIPLETS where spike times are the same, but trial IDs are shuffled

shuffletimes_T = 100;
shuffleLength_T = zeros(1,length(tripletSave));
shuffleLength_T_1_2_no3 = zeros(1,length(tripletSave));
for h=1:shuffletimes_T
    tic
    for i=1:size(tripletSave,1)
        
        firstcell=cell2mat(tripletSave(i,1));
        secondcell=cell2mat(tripletSave(i,2));
        thirdcell=cell2mat(tripletSave(i,3));
        
        newInd1L = lefttrials_choice(randperm(length(lefttrials_choice)));
        newInd1R = righttrials_choice(randperm(length(righttrials_choice)));
        newInd2L = lefttrials_choice(randperm(length(lefttrials_choice)));
        newInd2R = righttrials_choice(randperm(length(righttrials_choice)));
        newInd3L = lefttrials_choice(randperm(length(lefttrials_choice)));
        newInd3R = righttrials_choice(randperm(length(righttrials_choice)));
        
        savedShuffle_T = [];
        savedShuffle_T_1_2_no3 = [];
        
        Lc = 1;
        Rc = 1;
        for j=1:length(alltrials_choice)
            
            
            if alltrials_choice(j)==1
                trial1 = saveAll_preprocess(newInd1L(Lc)).digitized_eventStart(:,firstcell);
                trial2 = saveAll_preprocess(newInd2L(Lc)).digitized_eventStart(:,secondcell);
                trial3 = saveAll_preprocess(newInd3L(Lc)).digitized_eventStart(:,thirdcell);
                Lc=Lc+1;
            else
                trial1 = saveAll_preprocess(newInd1R(Rc)).digitized_eventStart(:,firstcell);
                trial2 = saveAll_preprocess(newInd2R(Rc)).digitized_eventStart(:,secondcell);
                trial3 = saveAll_preprocess(newInd3R(Rc)).digitized_eventStart(:,thirdcell);
                Rc=Rc+1;
            end
            
            cell1times = find(trial1==1);
            cell2times = find(trial2==1);
            cell3times = find(trial3==1);
            
            % identify the triplets
            if ~isempty(cell1times) && ~isempty(cell2times) && ~isempty(cell3times)
                if sum(cell1times(1)<cell2times)>0 && sum(cell3times(end)>cell2times)>0
                    savedShuffle_T = [savedShuffle_T alltrials_choice(j)];
                    
                    % case where 3 isn't after 2
                elseif sum(cell1times(1)<cell2times)>0 && sum(cell3times(end)>cell2times)==0
                    savedShuffle_T_1_2_no3 = [savedShuffle_T_1_2_no3 alltrials_choice(j)];
                end
                
                % case where 3 is empty but 2 after 1
            elseif ~isempty(cell1times) && ~isempty(cell2times) && isempty(cell3times)
                if sum(cell1times(1)<cell2times)>0
                    savedShuffle_T_1_2_no3 = [savedShuffle_T_1_2_no3 alltrials_choice(j)];
                end
            end
            
        end
        
        
        savedShuffles_T{h,i} = savedShuffle_T;
        savedShuffles_T_1_2_no3{h,i} = savedShuffle_T_1_2_no3;
        shuffleLength_T(1,i) = shuffleLength_T(1,i)+length(savedShuffle_T);
        shuffleLength_T_1_2_no3(1,i) = shuffleLength_T_1_2_no3(1,i)+length(savedShuffle_T_1_2_no3);
        
    end
    h
    toc
end

for i=1:size(savedShuffles_T,2)
    savedT3 = [];
    savedT3_1_2_no3 = [];
    for j=1:size(savedShuffles_T,1)
        savedT3 = [savedT3 savedShuffles_T{j,i}];
        savedT3_1_2_no3 = [savedT3_1_2_no3 savedShuffles_T_1_2_no3{j,i}];
    end
    savedShuffles2_T{i} = savedT3;
    savedShuffles2_T_1_2_no3{i} = savedT3_1_2_no3;
end

for i=1:size(savedShuffles2_T,2)
    
    savedShuffles3_T(i) = length(find(savedShuffles2_T{i}==1))/length(savedShuffles2_T{i});
    savedShuffles3_T_1_2_no3(i) = length(find(savedShuffles2_T_1_2_no3{i}==1))/length(savedShuffles2_T_1_2_no3{i});
    
end


% Code that calculated mean and std of length of shuffled data
savedShuffles_T_Length = cellfun(@length,savedShuffles_T);
means1 = mean(savedShuffles_T_Length);
stds1 = std(savedShuffles_T_Length);



savedShuffles3_T1 = num2cell(savedShuffles3_T);
[saveAll_triplets.triplet_shuffle] = savedShuffles3_T1{:};

savedShuffles3_T_1_2_no31 = num2cell(savedShuffles3_T_1_2_no3);
[saveAll_triplets.trip1_2_no3_shuffle] = savedShuffles3_T_1_2_no31{:};

savedSig1 = num2cell(sig1T(1,:));
[saveAll_triplets.sig_predict] = savedSig1{:};

means2 = num2cell(means1);
[saveAll_triplets.shuffle_length_mean] = means2{:};

stds2 = num2cell(stds1);
[saveAll_triplets.shuffle_length_std] = stds2{:};

for i=1:length(saveAll_triplets)
    
    if saveAll_triplets(i).prediction<.5
        saveAll_triplets(i).direction=1;
    else
        saveAll_triplets(i).direction=0;
    end
    
end

%% Calculate similarity of trials in which doublets show up versus random selection of the same # of L and R trials
% changed the way that seqanal_trialdistance works, so this doesn't work
% anymore

% 
% trialScores = seqanal_trialscore(score_subtrials, linspace(0,200,5));
% 
% for i=1:length(saveAll_triplets)
%     
%     exp_trials = saveAll_triplets(i).trials_appear;
%     exp_distance = seqanal_trialdistance(exp_trials, trialScores);
%     
%     num_left = length(find(saveAll_triplets(i).trials_choices==1));
%     num_right = length(find(saveAll_triplets(i).trials_choices==2));
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
% [saveAll_triplets.trials_distance] = exp_scoreCell{:};
% rand_score_meanCell = num2cell(rand_score_mean);
% [saveAll_triplets.random_distance_mean] = rand_score_meanCell{:};
% rand_score_stdCell = num2cell(rand_score_std);
% [saveAll_triplets.random_distance_std] = rand_score_stdCell{:};
% score_sigCell = num2cell(score_sig);
% [saveAll_triplets.score_sig] = score_sigCell{:};


%% Start Sig calculation


saveAll_triplets_Sig = saveAll_triplets;
saveAll_triplets_Sig(find([saveAll_triplets(:).sig_predict]==0))=[];
saveAll_triplets_Sig_Left = saveAll_triplets_Sig;
saveAll_triplets_Sig_Left(find([saveAll_triplets_Sig_Left(:).direction]==1))=[];
saveAll_triplets_Sig_Right = saveAll_triplets_Sig;
saveAll_triplets_Sig_Right(find([saveAll_triplets_Sig_Right(:).direction]==0))=[];

saveAll_triplets_great4 = saveAll_triplets;
saveAll_triplets_great4(find([saveAll_triplets(:).number_triplets]<4))=[];


saveAll_triplets_NotSig = saveAll_triplets;
saveAll_triplets_NotSig(find([saveAll_triplets(:).sig_predict]==1))=[];
saveAll_triplets_NotSig(find([saveAll_triplets_NotSig(:).number_triplets]<4))=[];


outputTriplets.saveAll_triplets = saveAll_triplets;
outputTriplets.saveAll_triplets_Sig = saveAll_triplets_Sig;
outputTriplets.saveAll_triplets_great4 = saveAll_triplets_great4;
outputTriplets.saveAll_triplets_Sig_Left = saveAll_triplets_Sig_Left;
outputTriplets.saveAll_triplets_Sig_Right = saveAll_triplets_Sig_Right;
outputTriplets.saveAll_triplets_NotSig = saveAll_triplets_NotSig;

%% Do the plotting

findTriplets_20190509_plotter(outputTriplets);
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(3,3,1)
% plot([saveAll_triplets_Sig_Left.triplet_shuffle],'LineWidth',2);
% hold on;
% plot([saveAll_triplets_Sig_Left.prediction],'LineWidth',2);
% legend('Scamble','Real');
% title('Significant Left-Preferring Triplets');
% xlabel('Triplet');
% ylabel('Fraction Trials Going Left');
% 
% subplot(3,3,2)
% plot([saveAll_triplets_Sig_Right.triplet_shuffle],'LineWidth',2);
% hold on;
% plot([saveAll_triplets_Sig_Right.prediction],'LineWidth',2);
% legend('Scamble','Real');
% title('Significant Right-Preferring Triplets');
% xlabel('Triplet');
% ylabel('Fraction Trials Going Left');
% 
% subplot(3,3,3)
% plot([saveAll_triplets_great5.shuffle_length_mean],'LineWidth',2);
% hold on;
% plot([saveAll_triplets_great5.number_triplets],'LineWidth',2);
% ylim([0 90]);
% legend('Scamble','Real');
% xlabel('Triplet (>5 Instances)');
% ylabel('# Trials with Triplet');
% 
% subplot(3,3,4)
% ht1 = histogram([saveAll_triplets_Sig.prediction]);
% ht1.BinWidth = 0.05;
% hold on;
% ht2 = histogram([saveAll_triplets_NotSig.prediction]);
% ht2.BinWidth = 0.05;
% xlabel('Fraction Trials Going Left');
% title('Significant Triplets - Histogram');
% ylabel('# of Triplets');
% 
% subplot(3,3,5)
% histogram(sort([saveAll_triplets_Sig_Left.prediction]-[saveAll_triplets_Sig_Left.triplet_shuffle]),[-1:.05:1]);
% hold on
% histogram(sort([saveAll_triplets_Sig_Right.prediction]-[saveAll_triplets_Sig_Right.triplet_shuffle]),[-1:.05:1]);
% xlim([-.8 .8]);
% legend('Left Doublets', 'Right Doublets');
% title('Real - Shuffle Prediction');
% xlabel('Doublet');
% ylabel('Difference Fraction Trials Going Left');
% 
% subplot(3,3,6)
% histogram([saveAll_triplets_Sig_Right.cell1and2only_predict], [0:.05:1])
% hold on
% histogram([saveAll_triplets_Sig_Right.prediction], [0:.05:1])
% histogram([saveAll_triplets_Sig_Right.cell3only_predict],[0:.05:1])
% legend('1 and 2 only','1, 2, and 3','3 only');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% title('Triplets - Right Predicting');
% 
% 
% subplot(3,3,9)
% histogram([saveAll_triplets_Sig_Left.cell1and2only_predict], [0:.05:1])
% hold on
% histogram([saveAll_triplets_Sig_Left.prediction], [0:.05:1])
% histogram([saveAll_triplets_Sig_Left.cell3only_predict],[0:.05:1])
% legend('1 and 2 only','1, 2, and 3','3 only');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% title('Triplets - Left Predicting');
% 
% subplot(3,3,7)
% histogram([saveAll_triplets_Sig_Right(:).cell1and2only_predict],10);
% hold on;
% histogram([saveAll_triplets_Sig_Right(:).cell2and3only_predict],10);
% % ylim([0 1500]);
% legend('Cell 1 and 2 of Triplet ONLY', 'Cell 2 and 3 of Triplet ONLY','Location','northwest');
% title('Triplets - Right Predicting');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% 
% subplot(3,3,8)
% histogram([saveAll_triplets_Sig_Left(:).cell1and2only_predict],10);
% hold on;
% histogram([saveAll_triplets_Sig_Left(:).cell2and3only_predict],10);
% % ylim([0 1500]);
% legend('Cell 1 and 2 of Triplet ONLY', 'Cell 2 and 3 of Triplet ONLY','Location','northwest');
% title('Triplets - Left Predicting');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
