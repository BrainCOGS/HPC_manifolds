function outputBehavior = extractBehavior(fname)

%% Function for extracting behavioral variables from the trial structure

% Used to generate the structure to input to functions to calculate the
% psychometric curve and logistic regression

% fname is the modeling.mat file that contains score and trial structures

% **********************************************************************
% Use CatStructFields to concatenate the outputBehavior of multiple animals
% **********************************************************************


%% Get mouseID and date from the filename
[mouseID, date, isModeling] = parseModelingName(fname);

%% get main trials from modeling file
if isModeling==1
    load(fname);
    mainTrials = trial([score.trial.mainTrial]==1);
    
%% get main trials from log file
else
    load(fname);
    
    for i=1:length(log.block)
        tempNums = num2cell([log.block(i).firstTrial:(log.block(i).firstTrial-1) + (length(log.block(i).trial))]);
        [log.block(i).trial.globalTrialID] = tempNums{:};
    end
    
    mainBlocks = find([log.block.mazeID]==11 | [log.block.mazeID]==7);
    mainTrials = [log.block(mainBlocks).trial];
    mainTrials = [mainTrials.globalTrialID];
    allTrials = [log.block(:).trial];
end

%%

for i=1:length(allTrials)
   tempPos_L{i} = allTrials(i).cuePos{1};
   tempPos_R{i} = allTrials(i).cuePos{2};
   numL = sum(allTrials(i).cueCombo(1,:));
   numR = sum(allTrials(i).cueCombo(2,:));
   tempRminusL(i) = numR-numL;
   tempRplusL(i) = numR+numL;
   
   position = allTrials(i).position;
   yposition = position(:,2);
   collision = allTrials(i).collision;
   collision1 = find(collision==1);
   if sum(collision1)==0
       ypos_collide(i) = 300;
   else
       ypos_collide(i) = yposition(collision1(1));
   end
   
   % Right choice is 1 and left choice is 0
   tempChoice = allTrials(i).choice;
   if tempChoice==1
       tempChoice2(i) = 0;
   elseif tempChoice==2
       tempChoice2(i) = 1;
   else
       tempChoice2(i) = NaN;
       % disp(['Choice was NaN in: ' fname]);
   end
   
   tempType = allTrials(i).trialType;
   if tempType==1
       tempType2(i) = 0;
   elseif tempType==2
       tempType2(i) = 1;
   else
       tempType2(i) = NaN;
   end
   
   if tempType2(i)==tempChoice2(i)
       tempChoiceCorrect(i) = 1;
   else
       tempChoiceCorrect(i) = 0;
   end
   
   if i==1
       tempPriorChoice(i) = NaN;
       tempPriorCorrect(i) = NaN;
   else
       tempPriorChoice(i) = tempChoice2(i-1);
       tempPriorCorrect(i) = tempChoiceCorrect(i-1);
   end
   
end

% outputBehavior.mouseID = repmat(mouseID, 1, length(mainTrials));
% outputBehavior.date = repmat(date, 1, length(mainTrials));
% outputBehavior.sessionID = repmat(1, 1, length(mainTrials));

%% Get rid of trials where the choice was not made
% keepData = ~isnan(tempChoice2);
keepData = ~isnan(tempChoice2) & ~isnan(tempPriorChoice) & ~isnan(tempPriorCorrect) & ismember([allTrials.globalTrialID], mainTrials);

outputBehavior.mouseID = repmat(mouseID, 1, sum(keepData==1));
outputBehavior.date = repmat(date, 1, sum(keepData==1));
outputBehavior.sessionID = repmat(1, 1, sum(keepData==1));
outputBehavior.cuePos_L = tempPos_L(keepData);
outputBehavior.cuePos_R = tempPos_R(keepData);
outputBehavior.choice   = tempChoice2(keepData);
outputBehavior.trialType = tempType2(keepData);
outputBehavior.choiceCorrect = tempChoiceCorrect(keepData);
outputBehavior.priorChoice = tempPriorChoice(keepData);
outputBehavior.priorCorrect = tempPriorCorrect(keepData);
outputBehavior.nCues_RminusL = tempRminusL(keepData);
outputBehavior.nCues_RplusL = tempRplusL(keepData);
outputBehavior.ypos_collide = ypos_collide(keepData);



%% Log of fnames and commands
% 
% fname_E74_20181205 = '\\bucket.pni.princeton.edu\Bezos-center\RigData\scope\bay3\edward\PoissonTowers_5\imaging\E74\20181205\E74_20181205_40per_userSetSD5minDur0.modeling.mat';
% fname_E74_20181206 = '\\bucket.pni.princeton.edu\Bezos-center\RigData\scope\bay3\edward\PoissonTowers_5\imaging\E74\20181206\E74_20181206_40per_userSetSD5minDur0.modelingCLEAN.mat';
% fname_E74_20181207 = '\\bucket.pni.princeton.edu\Bezos-center\RigData\scope\bay3\edward\PoissonTowers_5\imaging\E74\20181207\E74_20181207_50per_userSetSD5minDur0.modeling.mat';
% fname_E74_20181210 = '\\bucket.pni.princeton.edu\Bezos-center\RigData\scope\bay3\edward\PoissonTowers_5\imaging\E74\20181210\E74_20181210_50per_userSetSD5minDur0.modeling.mat';
% fname_E74_20181211 = '\\bucket.pni.princeton.edu\Bezos-center\RigData\scope\bay3\edward\PoissonTowers_5\imaging\E74\20181211\E74_20181211_50per_userSetSD5minDur0.modeling.mat';
% fname_E74_20181212 = '\\bucket.pni.princeton.edu\Bezos-center\RigData\scope\bay3\edward\PoissonTowers_5\imaging\E74\20181212\E74_20181212_50per_userSetSD5minDur0.modeling.mat';
% outputBehavior_E74_20181205 = extractBehavior(fname_E74_20181205);
% outputBehavior_E74_20181206 = extractBehavior(fname_E74_20181206);
% outputBehavior_E74_20181207 = extractBehavior(fname_E74_20181207);
% outputBehavior_E74_20181210 = extractBehavior(fname_E74_20181210);
% outputBehavior_E74_20181211 = extractBehavior(fname_E74_20181211);
% outputBehavior_E74_20181212 = extractBehavior(fname_E74_20181212);
% 
% 
% 
% 
% fname_E22_20170215 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170215.mat';
% fname_E22_20170216 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170216.mat';
% fname_E22_20170217 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170217.mat';
% fname_E22_20170221 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170221.mat';
% fname_E22_20170222 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170222.mat';
% fname_E22_20170223 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170223.mat';
% fname_E22_20170224 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170224.mat';
% fname_E22_20170227 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170227.mat';
% fname_E22_20170228 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170228.mat';
% fname_E22_20170308 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170308.mat';
% fname_E22_20170309 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170309.mat';
% fname_E22_20170320 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170320.mat';
% fname_E22_20170321 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170321.mat';
% fname_E22_20170403 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170403.mat';
% fname_E22_20170412 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170412.mat';
% fname_E22_20170424 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170424.mat';
% 
% outputBehavior_E22_20170215 = extractBehavior(fname_E22_20170215);
% outputBehavior_E22_20170216 = extractBehavior(fname_E22_20170216);
% outputBehavior_E22_20170217 = extractBehavior(fname_E22_20170217);
% outputBehavior_E22_20170221 = extractBehavior(fname_E22_20170221);
% outputBehavior_E22_20170222 = extractBehavior(fname_E22_20170222);
% outputBehavior_E22_20170223 = extractBehavior(fname_E22_20170223);
% outputBehavior_E22_20170224 = extractBehavior(fname_E22_20170224);
% outputBehavior_E22_20170227 = extractBehavior(fname_E22_20170227);
% outputBehavior_E22_20170228 = extractBehavior(fname_E22_20170228);
% outputBehavior_E22_20170308 = extractBehavior(fname_E22_20170308);
% outputBehavior_E22_20170309 = extractBehavior(fname_E22_20170309);
% outputBehavior_E22_20170320 = extractBehavior(fname_E22_20170320);
% outputBehavior_E22_20170321 = extractBehavior(fname_E22_20170321);
% outputBehavior_E22_20170403 = extractBehavior(fname_E22_20170403);
% outputBehavior_E22_20170412 = extractBehavior(fname_E22_20170412);
% outputBehavior_E22_20170424 = extractBehavior(fname_E22_20170424);
% 
% 
% 
% fname_E39_20171011 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171011.mat';
% fname_E39_20171012 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171012.mat';
% fname_E39_20171102 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171102.mat';
% fname_E39_20171103 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171103.mat';
% fname_E39_20171110 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171110.mat';
% 
% outputBehavior_E39_20171011 = extractBehavior(fname_E39_20171011);
% outputBehavior_E39_20171012 = extractBehavior(fname_E39_20171012);
% outputBehavior_E39_20171102 = extractBehavior(fname_E39_20171102);
% outputBehavior_E39_20171103 = extractBehavior(fname_E39_20171103);
% outputBehavior_E39_20171110 = extractBehavior(fname_E39_20171110);
% 
% 
% 
% 
% fname_E43_20170720 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170720.mat';
% fname_E43_20170721 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170721.mat';
% fname_E43_20170724 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170724.mat';
% fname_E43_20170726 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170726.mat';
% fname_E43_20170801 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170801.mat';
% fname_E43_20170802 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170802.mat';
% fname_E43_20170803 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170803.mat';
% fname_E43_20170804 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170804.mat';
% fname_E43_20170807 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170807.mat';
% fname_E43_20170811 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170811.mat';
% fname_E43_20170815 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170815.mat';
% fname_E43_20170816 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170816.mat';
% fname_E43_20170817 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170817.mat';
% fname_E43_20170818 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170818.mat';
% fname_E43_20170821 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170821.mat';
% fname_E43_20170825 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170825.mat';
% fname_E43_20170829 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170829.mat';
% fname_E43_20170831 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170831.mat';
% fname_E43_20170905 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170905.mat';
% fname_E43_20170911 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170911.mat';
% fname_E43_20170912 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170912.mat';
% 
% outputBehavior_E43_20170720 = extractBehavior(fname_E43_20170720);
% outputBehavior_E43_20170721 = extractBehavior(fname_E43_20170721);
% outputBehavior_E43_20170724 = extractBehavior(fname_E43_20170724);
% outputBehavior_E43_20170726 = extractBehavior(fname_E43_20170726);
% outputBehavior_E43_20170801 = extractBehavior(fname_E43_20170801);
% outputBehavior_E43_20170802 = extractBehavior(fname_E43_20170802);
% outputBehavior_E43_20170803 = extractBehavior(fname_E43_20170803);
% outputBehavior_E43_20170804 = extractBehavior(fname_E43_20170804);
% outputBehavior_E43_20170807 = extractBehavior(fname_E43_20170807);
% outputBehavior_E43_20170811 = extractBehavior(fname_E43_20170811);
% outputBehavior_E43_20170815 = extractBehavior(fname_E43_20170815);
% outputBehavior_E43_20170816 = extractBehavior(fname_E43_20170816);
% outputBehavior_E43_20170817 = extractBehavior(fname_E43_20170817);
% outputBehavior_E43_20170818 = extractBehavior(fname_E43_20170818);
% outputBehavior_E43_20170821 = extractBehavior(fname_E43_20170821);
% outputBehavior_E43_20170825 = extractBehavior(fname_E43_20170825);
% outputBehavior_E43_20170829 = extractBehavior(fname_E43_20170829);
% outputBehavior_E43_20170831 = extractBehavior(fname_E43_20170831);
% outputBehavior_E43_20170905 = extractBehavior(fname_E43_20170905);
% outputBehavior_E43_20170911 = extractBehavior(fname_E43_20170911);
% outputBehavior_E43_20170912 = extractBehavior(fname_E43_20170912);
% 
% 
% fname_E44_20170914 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170914.mat';
% fname_E44_20170915 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170915.mat';
% fname_E44_20170918 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170918.mat';
% fname_E44_20170919 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170919.mat';
% fname_E44_20170920 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170920.mat';
% fname_E44_20170921 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170921.mat';
% fname_E44_20170922 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170922.mat';
% fname_E44_20170925 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170925.mat';
% fname_E44_20171010 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171010.mat';
% fname_E44_20171011 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171011.mat';
% fname_E44_20171012 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171012.mat';
% fname_E44_20171013 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171013.mat';
% fname_E44_20171017 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171017.mat';
% fname_E44_20171018 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171018.mat';
% fname_E44_20171019 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171019.mat';
% fname_E44_20171020 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171020.mat';
% fname_E44_20171024 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171024.mat';
% fname_E44_20171025 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171025.mat';
% fname_E44_20171108 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171108.mat';
% 
% outputBehavior_E44_20170914 = extractBehavior(fname_E44_20170914);
% outputBehavior_E44_20170915 = extractBehavior(fname_E44_20170915);
% outputBehavior_E44_20170918 = extractBehavior(fname_E44_20170918);
% outputBehavior_E44_20170919 = extractBehavior(fname_E44_20170919);
% outputBehavior_E44_20170920 = extractBehavior(fname_E44_20170920);
% outputBehavior_E44_20170921 = extractBehavior(fname_E44_20170921);
% outputBehavior_E44_20170922 = extractBehavior(fname_E44_20170922);
% outputBehavior_E44_20170925 = extractBehavior(fname_E44_20170925);
% outputBehavior_E44_20171010 = extractBehavior(fname_E44_20171010);
% outputBehavior_E44_20171011 = extractBehavior(fname_E44_20171011);
% outputBehavior_E44_20171012 = extractBehavior(fname_E44_20171012);
% outputBehavior_E44_20171013 = extractBehavior(fname_E44_20171013);
% outputBehavior_E44_20171017 = extractBehavior(fname_E44_20171017);
% outputBehavior_E44_20171018 = extractBehavior(fname_E44_20171018);
% outputBehavior_E44_20171019 = extractBehavior(fname_E44_20171019);
% outputBehavior_E44_20171020 = extractBehavior(fname_E44_20171020);
% outputBehavior_E44_20171024 = extractBehavior(fname_E44_20171024);
% outputBehavior_E44_20171025 = extractBehavior(fname_E44_20171025);
% outputBehavior_E44_20171108 = extractBehavior(fname_E44_20171108);
% 
% 
% 
% 
% fname_E47_20170927 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20170927.mat';
% fname_E47_20171003 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171003.mat';
% fname_E47_20171004 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171004.mat';
% fname_E47_20171005 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171005.mat';
% fname_E47_20171006 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171006.mat';
% fname_E47_20171010 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171010.mat';
% fname_E47_20171011 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171011.mat';
% fname_E47_20171012 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171012.mat';
% fname_E47_20171013 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171013.mat';
% fname_E47_20171018 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171018.mat';
% outputBehavior_E47_20170927 = extractBehavior(fname_E47_20170927);
% outputBehavior_E47_20171003 = extractBehavior(fname_E47_20171003);
% outputBehavior_E47_20171004 = extractBehavior(fname_E47_20171004);
% outputBehavior_E47_20171005 = extractBehavior(fname_E47_20171005);
% outputBehavior_E47_20171006 = extractBehavior(fname_E47_20171006);
% outputBehavior_E47_20171010 = extractBehavior(fname_E47_20171010);
% outputBehavior_E47_20171011 = extractBehavior(fname_E47_20171011);
% outputBehavior_E47_20171012 = extractBehavior(fname_E47_20171012);
% outputBehavior_E47_20171013 = extractBehavior(fname_E47_20171013);
% outputBehavior_E47_20171018 = extractBehavior(fname_E47_20171018);
% 
% 
% 
% fname_E48_20170807 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170807.mat';
% fname_E48_20170808 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170808.mat';
% fname_E48_20170809 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170809.mat';
% fname_E48_20170810 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170810.mat';
% fname_E48_20170811 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170811.mat';
% fname_E48_20170814 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170814.mat';
% fname_E48_20170815 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170815.mat';
% fname_E48_20170816 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170816.mat';
% fname_E48_20170817 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170817.mat';
% fname_E48_20170823 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170823.mat';
% fname_E48_20170824 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170824.mat';
% fname_E48_20170828 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170828.mat';
% fname_E48_20170830 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170830.mat';
% fname_E48_20170831 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170831.mat';
% fname_E48_20170905 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170905.mat';
% fname_E48_20170906 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170906.mat';
% fname_E48_20170907 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170907.mat';
% fname_E48_20170908 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170908.mat';
% fname_E48_20170911 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170911.mat';
% fname_E48_20170914 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170914.mat';
% fname_E48_20170915 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170915.mat';
% fname_E48_20170918 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170918.mat';
% fname_E48_20170919 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170919.mat';
% fname_E48_20170920 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170920.mat';
% fname_E48_20171004 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20171004.mat';
% 
% outputBehavior_E48_20170807 = extractBehavior(fname_E48_20170807);
% outputBehavior_E48_20170808 = extractBehavior(fname_E48_20170808);
% outputBehavior_E48_20170809 = extractBehavior(fname_E48_20170809);
% outputBehavior_E48_20170810 = extractBehavior(fname_E48_20170810);
% outputBehavior_E48_20170811 = extractBehavior(fname_E48_20170811);
% outputBehavior_E48_20170814 = extractBehavior(fname_E48_20170814);
% outputBehavior_E48_20170815 = extractBehavior(fname_E48_20170815);
% outputBehavior_E48_20170816 = extractBehavior(fname_E48_20170816);
% outputBehavior_E48_20170817 = extractBehavior(fname_E48_20170817);
% outputBehavior_E48_20170823 = extractBehavior(fname_E48_20170823);
% outputBehavior_E48_20170824 = extractBehavior(fname_E48_20170824);
% outputBehavior_E48_20170828 = extractBehavior(fname_E48_20170828);
% outputBehavior_E48_20170830 = extractBehavior(fname_E48_20170830);
% outputBehavior_E48_20170831 = extractBehavior(fname_E48_20170831);
% outputBehavior_E48_20170905 = extractBehavior(fname_E48_20170905);
% outputBehavior_E48_20170906 = extractBehavior(fname_E48_20170906);
% outputBehavior_E48_20170907 = extractBehavior(fname_E48_20170907);
% outputBehavior_E48_20170908 = extractBehavior(fname_E48_20170908);
% outputBehavior_E48_20170911 = extractBehavior(fname_E48_20170911);
% outputBehavior_E48_20170914 = extractBehavior(fname_E48_20170914);
% outputBehavior_E48_20170915 = extractBehavior(fname_E48_20170915);
% outputBehavior_E48_20170918 = extractBehavior(fname_E48_20170918);
% outputBehavior_E48_20170919 = extractBehavior(fname_E48_20170919);
% outputBehavior_E48_20170920 = extractBehavior(fname_E48_20170920);
% outputBehavior_E48_20171004 = extractBehavior(fname_E48_20171004);
% 
% 
% 
% fname_E52_20171221 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20171221.mat';
% fname_E52_20180103 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180103.mat';
% fname_E52_20180105 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180105.mat';
% fname_E52_20180108 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180108.mat';
% fname_E52_20180116 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180116.mat';
% fname_E52_20180117 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180117.mat';
% fname_E52_20180118 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180118.mat';
% fname_E52_20180126 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180126.mat';
% fname_E52_20180127 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180127.mat';
% fname_E52_20180129 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180129.mat';
% fname_E52_20180130 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180130.mat';
% fname_E52_20180201 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180201.mat';
% fname_E52_20180202 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180202.mat';
% fname_E52_20180205 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E52\PoissonBlocksReboot3_cohort3_Bezos3_E52_T_20180205.mat';
% 
% outputBehavior_E52_20171221 = extractBehavior(fname_E52_20171221);
% outputBehavior_E52_20180103 = extractBehavior(fname_E52_20180103);
% outputBehavior_E52_20180105 = extractBehavior(fname_E52_20180105);
% outputBehavior_E52_20180108 = extractBehavior(fname_E52_20180108);
% outputBehavior_E52_20180116 = extractBehavior(fname_E52_20180116);
% outputBehavior_E52_20180117 = extractBehavior(fname_E52_20180117);
% outputBehavior_E52_20180118 = extractBehavior(fname_E52_20180118);
% outputBehavior_E52_20180126 = extractBehavior(fname_E52_20180126);
% outputBehavior_E52_20180127 = extractBehavior(fname_E52_20180127);
% outputBehavior_E52_20180129 = extractBehavior(fname_E52_20180129);
% outputBehavior_E52_20180130 = extractBehavior(fname_E52_20180130);
% outputBehavior_E52_20180201 = extractBehavior(fname_E52_20180201);
% outputBehavior_E52_20180202 = extractBehavior(fname_E52_20180202);
% outputBehavior_E52_20180205 = extractBehavior(fname_E52_20180205);
% 
% 
% fname_E62_20171220 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20171220.mat';
% fname_E62_20171221 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20171221.mat';
% fname_E62_20180108 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180108.mat';
% fname_E62_20180109 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180109.mat';
% fname_E62_20180111 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180111.mat';
% fname_E62_20180112 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180112.mat';
% fname_E62_20180115 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180115.mat';
% fname_E62_20180117 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180117.mat';
% fname_E62_20180118 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180118.mat';
% fname_E62_20180119 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180119.mat';
% fname_E62_20180122 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180122.mat';
% fname_E62_20180123 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180123.mat';
% fname_E62_20180124 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180124.mat';
% fname_E62_20180125 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180125.mat';
% fname_E62_20180126 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180126.mat';
% fname_E62_20180129 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180129.mat';
% fname_E62_20180130 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180130.mat';
% fname_E62_20180131 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E62\PoissonBlocksReboot4_cohort4_Bezos3_E62_T_20180131.mat';
% 
% outputBehavior_E62_20171220 = extractBehavior(fname_E62_20171220);
% outputBehavior_E62_20171221 = extractBehavior(fname_E62_20171221);
% outputBehavior_E62_20180108 = extractBehavior(fname_E62_20180108);
% outputBehavior_E62_20180109 = extractBehavior(fname_E62_20180109);
% outputBehavior_E62_20180111 = extractBehavior(fname_E62_20180111);
% outputBehavior_E62_20180112 = extractBehavior(fname_E62_20180112);
% outputBehavior_E62_20180115 = extractBehavior(fname_E62_20180115);
% outputBehavior_E62_20180117 = extractBehavior(fname_E62_20180117);
% outputBehavior_E62_20180118 = extractBehavior(fname_E62_20180118);
% outputBehavior_E62_20180119 = extractBehavior(fname_E62_20180119);
% outputBehavior_E62_20180122 = extractBehavior(fname_E62_20180122);
% outputBehavior_E62_20180123 = extractBehavior(fname_E62_20180123);
% outputBehavior_E62_20180124 = extractBehavior(fname_E62_20180124);
% outputBehavior_E62_20180125 = extractBehavior(fname_E62_20180125);
% outputBehavior_E62_20180126 = extractBehavior(fname_E62_20180126);
% outputBehavior_E62_20180129 = extractBehavior(fname_E62_20180129);
% outputBehavior_E62_20180130 = extractBehavior(fname_E62_20180130);
% outputBehavior_E62_20180131 = extractBehavior(fname_E62_20180131);
% 
% 
% fname_E65_20180131 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180131.mat';
% fname_E65_20180201 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180201.mat';
% fname_E65_20180202 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180202.mat';
% fname_E65_20180205 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180205.mat';
% fname_E65_20180206 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180206.mat';
% fname_E65_20180207 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180207.mat';
% fname_E65_20180208 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180208.mat';
% fname_E65_20180209 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180209.mat';
% fname_E65_20180212 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180212.mat';
% fname_E65_20180213 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180213.mat';
% fname_E65_20180214 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180214.mat';
% fname_E65_20180312 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180312.mat';
% 
% outputBehavior_E65_20180131 = extractBehavior(fname_E65_20180131);
% outputBehavior_E65_20180201 = extractBehavior(fname_E65_20180201);
% outputBehavior_E65_20180202 = extractBehavior(fname_E65_20180202);
% outputBehavior_E65_20180205 = extractBehavior(fname_E65_20180205);
% outputBehavior_E65_20180206 = extractBehavior(fname_E65_20180206);
% outputBehavior_E65_20180207 = extractBehavior(fname_E65_20180207);
% outputBehavior_E65_20180208 = extractBehavior(fname_E65_20180208);
% outputBehavior_E65_20180209 = extractBehavior(fname_E65_20180209);
% outputBehavior_E65_20180212 = extractBehavior(fname_E65_20180212);
% outputBehavior_E65_20180213 = extractBehavior(fname_E65_20180213);
% outputBehavior_E65_20180214 = extractBehavior(fname_E65_20180214);
% outputBehavior_E65_20180312 = extractBehavior(fname_E65_20180312);
% 
% 
% 
% 
% S = CatStructFields(2, outputBehavior_E22_20170215, ...
%                        outputBehavior_E22_20170216, ...
%                        outputBehavior_E22_20170217, ...
%                        outputBehavior_E22_20170221, ...
%                        outputBehavior_E22_20170222, ...
%                        outputBehavior_E22_20170223, ...
%                        outputBehavior_E22_20170224, ...
%                        outputBehavior_E22_20170227, ...
%                        outputBehavior_E22_20170228, ...
%                        outputBehavior_E22_20170308, ...
%                        outputBehavior_E22_20170309, ...
%                        outputBehavior_E22_20170320, ...
%                        outputBehavior_E22_20170321, ...
%                        outputBehavior_E22_20170403, ...
%                        outputBehavior_E22_20170412, ...
%                        outputBehavior_E39_20171011, ...
%                        outputBehavior_E39_20171012, ...
%                        outputBehavior_E39_20171102, ...
%                        outputBehavior_E39_20171103, ...
%                        outputBehavior_E39_20171110, ...
%                        outputBehavior_E43_20170720, ... 
%                        outputBehavior_E43_20170721, ... 
%                        outputBehavior_E43_20170724, ... 
%                        outputBehavior_E43_20170726, ... 
%                        outputBehavior_E43_20170801, ... 
%                        outputBehavior_E43_20170802, ... 
%                        outputBehavior_E43_20170803, ... 
%                        outputBehavior_E43_20170804, ... 
%                        outputBehavior_E43_20170807, ... 
%                        outputBehavior_E43_20170811, ... 
%                        outputBehavior_E43_20170815, ... 
%                        outputBehavior_E43_20170816, ... 
%                        outputBehavior_E43_20170818, ... 
%                        outputBehavior_E43_20170821, ... 
%                        outputBehavior_E43_20170825, ... 
%                        outputBehavior_E43_20170829, ... 
%                        outputBehavior_E43_20170831, ... 
%                        outputBehavior_E43_20170905, ... 
%                        outputBehavior_E43_20170911, ... 
%                        outputBehavior_E43_20170912, ... 
%                        outputBehavior_E44_20170914, ... 
%                        outputBehavior_E44_20170915, ... 
%                        outputBehavior_E44_20170918, ... 
%                        outputBehavior_E44_20170919, ... 
%                        outputBehavior_E44_20170920, ... 
%                        outputBehavior_E44_20170921, ... 
%                        outputBehavior_E44_20170922, ... 
%                        outputBehavior_E44_20170925, ... 
%                        outputBehavior_E44_20171010, ... 
%                        outputBehavior_E44_20171011, ... 
%                        outputBehavior_E44_20171012, ... 
%                        outputBehavior_E44_20171013, ... 
%                        outputBehavior_E44_20171017, ... 
%                        outputBehavior_E44_20171018, ... 
%                        outputBehavior_E44_20171019, ... 
%                        outputBehavior_E44_20171020, ... 
%                        outputBehavior_E44_20171024, ... 
%                        outputBehavior_E44_20171025, ... 
%                        outputBehavior_E44_20171108, ... 
%                        outputBehavior_E47_20170927, ...
%                        outputBehavior_E47_20171003, ...
%                        outputBehavior_E47_20171004, ...
%                        outputBehavior_E47_20171005, ...
%                        outputBehavior_E47_20171006, ...
%                        outputBehavior_E47_20171010, ...
%                        outputBehavior_E47_20171011, ...
%                        outputBehavior_E47_20171012, ...
%                        outputBehavior_E47_20171013, ...
%                        outputBehavior_E47_20171018, ...
%                        outputBehavior_E48_20170807, ...
%                        outputBehavior_E48_20170808, ...
%                        outputBehavior_E48_20170809, ...
%                        outputBehavior_E48_20170810, ...
%                        outputBehavior_E48_20170811, ...
%                        outputBehavior_E48_20170814, ...
%                        outputBehavior_E48_20170815, ...
%                        outputBehavior_E48_20170816, ...
%                        outputBehavior_E48_20170817, ...
%                        outputBehavior_E48_20170823, ...
%                        outputBehavior_E48_20170824, ...
%                        outputBehavior_E48_20170828, ...
%                        outputBehavior_E48_20170830, ...
%                        outputBehavior_E48_20170831, ...
%                        outputBehavior_E48_20170905, ...
%                        outputBehavior_E48_20170906, ...
%                        outputBehavior_E48_20170907, ...
%                        outputBehavior_E48_20170908, ...
%                        outputBehavior_E48_20170911, ...
%                        outputBehavior_E48_20170914, ...
%                        outputBehavior_E48_20170915, ...
%                        outputBehavior_E48_20170918, ...
%                        outputBehavior_E48_20170919, ...
%                        outputBehavior_E48_20170920, ...
%                        outputBehavior_E48_20171004, ...
%                        outputBehavior_E52_20171221, ...
%                        outputBehavior_E52_20180103, ...
%                        outputBehavior_E52_20180105, ...
%                        outputBehavior_E52_20180108, ...
%                        outputBehavior_E52_20180116, ...
%                        outputBehavior_E52_20180117, ...
%                        outputBehavior_E52_20180118, ...
%                        outputBehavior_E52_20180126, ...
%                        outputBehavior_E52_20180127, ...
%                        outputBehavior_E52_20180129, ...
%                        outputBehavior_E52_20180130, ...
%                        outputBehavior_E52_20180201, ...
%                        outputBehavior_E52_20180202, ...
%                        outputBehavior_E52_20180205, ...
%                        outputBehavior_E62_20171220, ...
%                        outputBehavior_E62_20171221, ...
%                        outputBehavior_E62_20180108, ...
%                        outputBehavior_E62_20180109, ...
%                        outputBehavior_E62_20180111, ...
%                        outputBehavior_E62_20180112, ...
%                        outputBehavior_E62_20180115, ...
%                        outputBehavior_E62_20180117, ...
%                        outputBehavior_E62_20180118, ...
%                        outputBehavior_E62_20180119, ...
%                        outputBehavior_E62_20180122, ...
%                        outputBehavior_E62_20180123, ...
%                        outputBehavior_E62_20180124, ...
%                        outputBehavior_E62_20180125, ...
%                        outputBehavior_E62_20180126, ...
%                        outputBehavior_E62_20180129, ...
%                        outputBehavior_E62_20180130, ...
%                        outputBehavior_E62_20180131, ...
%                        outputBehavior_E65_20180131, ...
%                        outputBehavior_E65_20180201, ...
%                        outputBehavior_E65_20180202, ...
%                        outputBehavior_E65_20180205, ...
%                        outputBehavior_E65_20180206, ...
%                        outputBehavior_E65_20180207, ...
%                        outputBehavior_E65_20180208, ...
%                        outputBehavior_E65_20180209, ...
%                        outputBehavior_E65_20180212, ...
%                        outputBehavior_E65_20180213, ...
%                        outputBehavior_E65_20180214, ...
%                        outputBehavior_E65_20180312, ...
%                        outputBehavior_E74_20181205, ...
%                        outputBehavior_E74_20181206, ...
%                        outputBehavior_E74_20181207, ...
%                        outputBehavior_E74_20181210, ...
%                        outputBehavior_E74_20181211, ...
%                        outputBehavior_E74_20181212 ...    
%                    );
% rc = logRegressionFromConcatLog2(S);

