function outputTriplets = findTriplets(outputDoublets)

config.rngInput = 1;
rng(config.rngInput);

saveAll_preprocess = outputDoublets.saveAll_preprocess;
saveAll_doublets   = outputDoublets.saveAll_doublets;
saveAll_basics     = outputDoublets.saveAll_basics;
score_subtrials    = outputDoublets.score_subtrials;
temp3a             = {outputDoublets.saveAll_preprocess.digitized_eventStart};


%% Find the triplets

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

% Get rid of triplets where first and third cell are the same
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


%% Get important triplet statistics

tic
for i=1:size(tripletSave,1)
    
    bothinorder1_T = cell2mat(tripletSave(i,5));
    bothinorder2_T = [alltrials_choice(bothinorder1_T)];
    saveOrdered_T(i) = length(find(bothinorder2_T==1))/length(bothinorder2_T);
    saveLength_T(i) = length(bothinorder2_T);
    saveChoices_T{i} = bothinorder2_T;
    
    onlycell1Choice = [alltrials_choice(cell2mat(tripletSave(i,6)))];
    onlycell2Choice = [alltrials_choice(cell2mat(tripletSave(i,7)))];
    onlycell3Choice = [alltrials_choice(cell2mat(tripletSave(i,8)))];
    
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

%% Create Structure for Triplets
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
    'cell1and3only_length',num2cell(saveLength_cell1and3)');


%% New shuffle FOR TRIPLETS where spike times are the same, but trial IDs are shuffled

shuffletimes_T = 100;
shuffleLength_T = zeros(1,length(tripletSave));
shuffleLength_T_1_2_no3 = zeros(1,length(tripletSave));

tic
for h=1:shuffletimes_T
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
end
toc

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


%% Set the direction of each doublet, 1 means right-predicting

for i=1:length(saveAll_triplets)
    
    if saveAll_triplets(i).prediction<.5
        saveAll_triplets(i).direction=1;
    else
        saveAll_triplets(i).direction=0;
    end
    
end

%% Start Sig calculation

saveAll_triplets_Sig = saveAll_triplets;
saveAll_triplets_Sig(find([saveAll_triplets(:).sig_predict]==0))=[];

saveAll_triplets_NotSig = saveAll_triplets;
saveAll_triplets_NotSig(find([saveAll_triplets(:).sig_predict]==1))=[];

saveAll_triplets_Sig_Left = saveAll_triplets_Sig;
saveAll_triplets_Sig_Left(find([saveAll_triplets_Sig_Left(:).direction]==1))=[];
saveAll_triplets_Sig_Right = saveAll_triplets_Sig;
saveAll_triplets_Sig_Right(find([saveAll_triplets_Sig_Right(:).direction]==0))=[];


%% Make the outputfile

outputTriplets.saveAll_triplets = saveAll_triplets;
outputTriplets.saveAll_triplets_Sig = saveAll_triplets_Sig;
outputTriplets.saveAll_triplets_Sig_Left = saveAll_triplets_Sig_Left;
outputTriplets.saveAll_triplets_Sig_Right = saveAll_triplets_Sig_Right;
outputTriplets.saveAll_triplets_NotSig = saveAll_triplets_NotSig;

