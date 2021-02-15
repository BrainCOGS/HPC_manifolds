function outputNonlinearDecoding = mind_nonlinearDecoding_dimX_All(fname, fname_mani, numFolds, regType, varType, taskType, shuffleToggle, manifoldToggle, roiList, dimEmbed)
% regType should be either GP or SVM or LLE
% varType should be variable that wants to be decoded
% taskType should be 'towers' or 'alternation', etc
% if shuffleToggle==1, shuffle the manifold data
% if manifoldToggle==1, use the manifold, if 0, use the dataDFF
% leave roiList empty

%% Set up config and argins
config.rng      =1;
rng(config.rng);

argins.fname = fname;
argins.fname_mani = fname_mani;
argins.numFolds   = numFolds;
argins.regType    = regType;
argins.varType    = varType;
argins.taskType   = taskType;
argins.shuffleToggle = shuffleToggle;
argins.manifoldToggle = manifoldToggle;
argins.roiList        = roiList;
argins.dimEmbed = dimEmbed;

outputNonlinearDecoding.config    = config;
outputNonlinearDecoding.argins    = argins;


%% Load the manifold

if ischar(fname_mani)
    load(fname_mani);
elseif isstring(fname_mani)
    load(fname_mani);
else
    outMind = fname_mani;
end

%% Use the correct extractVariables

if strcmp(taskType, 'alternation')==1 || strcmp(taskType, 'Alternation')==1
    nic_output = extractVariables('all', 2, 'goodTrials', 2, 0, 0, 5, [0 0], fname,'none','alternation', 1, 1);
    
elseif strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [0 0], fname,'none', 'towers', 1, 1);
    
elseif strcmp(taskType, 'alternationJeff')==1 || strcmp(taskType, 'AlternationJeff')==1
    nic_output = extractVariables('all', 6, 'goodTrials', 2, 0, 0, 5, [0 0], fname,'none','alternation', 1, 1);
    
end


%% Take roiList if it exists

if ~isempty(roiList)
    nic_output.ROIactivities = nic_output.ROIactivities(:,roiList);
end

%% Get variables
behavioralVariables = nic_output.behavioralVariables;
ROIactivities = nic_output.ROIactivities;
ROIactivities_full = ROIactivities;
% if outMind.argins.downSample~=0
%     behavioralVariables = behavioralVariables(1:outMind.argins.downSample:end,:);
%     ROIactivities = ROIactivities(1:outMind.argins.downSample:end,:);
% end
behavioralVariables = behavioralVariables(outMind.Datarange,:);
ROIactivities = ROIactivities(outMind.Datarange,:);

trialn = behavioralVariables.Trial;
utrial = unique(trialn);


%% Manifold or ROI activity

if manifoldToggle==1
    
    % if outMind.argins.downSample>1
    %     sample = ROIactivities_full(outMind.config.Datarange,outMind.Neurons);
    % else
    sample = outMind.dataDFF;
    %end
    pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
    manifold3d = outMind.dat.allembed([outMind.dat.allembed.d]==dimEmbed).f2m.map.transform(pca_coords);
    
    dataUsed = manifold3d;
else
    dataUsed = ROIactivities;
    warndlg('Make sure that manifold loaded is not downsampled, otherwise it will use subset of data');
end

if shuffleToggle==1
    jeff_constant = 10;
    dataUsed = jeff_randomizer(dataUsed, jeff_constant);
end

%% Generate Crossvalidation
CV = generateCrossValSet_v2(behavioralVariables, numFolds);

dataUsed_full = dataUsed;

for i=1:numFolds
        
    if isnumeric(varType)
        roiout = ROIactivities(:,varType);
        trainResp1 = roiout(CV(i).trainLocations);
        testResp1 = roiout(CV(i).testLocations);
    elseif  strcmp(varType,'Evidence')
        trainResp1 = CV(i).trainData.Evidence;
        testResp1  = CV(i).testData.Evidence;
        varLetter  = 'E';
    elseif strcmp(varType, 'EvidenceSmooth')
        trainResp1 = CV(i).trainData.EvidenceSmooth;
        testResp1  = CV(i).testData.EvidenceSmooth;
        varLetter  = 'ES';
    elseif strcmp(varType, 'Position')
        trainResp1 = CV(i).trainData.Position;
        testResp1  = CV(i).testData.Position;
        varLetter  = 'Y';
    elseif strcmp(varType, 'ViewAngle')
        trainResp1 = CV(i).trainData.ViewAngle;
        testResp1  = CV(i).testData.ViewAngle;
    elseif strcmp(varType, 'Choice')
        trainResp1 = CV(i).trainData.Choice;
        testResp1  = CV(i).testData.Choice;
        posResp1  = CV(i).testData.Position;
    elseif strcmp(varType, 'ChoiceCorrect')
        trainResp1 = CV(i).trainData.ChoiceCorrect;
        testResp1  = CV(i).testData.ChoiceCorrect;
        posResp1  = CV(i).testData.Position;
    elseif strcmp(varType, 'PriorChoice')
        trainResp1 = CV(i).trainData.PriorChoice;
        testResp1  = CV(i).testData.PriorChoice;
        posResp1  = CV(i).testData.Position;
    elseif strcmp(varType, 'PriorCorrect')
        trainResp1 = CV(i).trainData.PriorCorrect;
        testResp1  = CV(i).testData.PriorCorrect;
        posResp1  = CV(i).testData.Position;
    elseif strcmp(varType, 'EvidenceR')
        trainResp1 = CV(i).trainData.EvidenceR;
        testResp1  = CV(i).testData.EvidenceR;
    elseif strcmp(varType, 'EvidenceL')
        trainResp1 = CV(i).trainData.EvidenceL;
        testResp1  = CV(i).testData.EvidenceL;
    elseif strcmp(varType, 'NumberOfCues')
        trainResp1 = CV(i).trainData.NumberOfCues;
        testResp1  = CV(i).testData.NumberOfCues;
    elseif strcmp(varType, 'Velocity')
        trainResp1 = CV(i).trainData.Velocity;
        testResp1  = CV(i).testData.Velocity;
    elseif strcmp(varType, 'EvidenceR-EvidenceL')
        trainResp1 = CV(i).trainData.EvidenceR;
        testResp1  = CV(i).testData.Evidence;
        trainResp2 = CV(i).trainData.EvidenceL;
    elseif strcmp(varType, 'PrevEvidence')
        trainResp1 = CV(i).trainData.PrevEvidence;
        testResp1  = CV(i).testData.PrevEvidence;
    end
    
    % Get the top 10% of Position cells if using Neural activity and no
    % roiList was provided
    if manifoldToggle~=1 && isempty(roiList)
        tempname = split(fname_mani,'\');
        tempname = tempname{end};
        tempname = tempname(1:3);
        if strcmp(varType,'Position')
            out_Data = getSkaggs([5 2], 'noLog', 1, CV(i).train', 2, [tempname '_' varLetter '_[5_2]_CV' num2str(i) '_%d.mat'], 1, fname, 2, {'Position'}, {[0:10:300]}, {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
        elseif strcmp(varType,'Evidence')
            out_Data = getSkaggs([5 2], 'noLog', 1, CV(i).train', 2, [tempname '_' varLetter '_[5_2]_CV' num2str(i) '_%d.mat'], 1, fname, 2, {'Evidence'}, {[-15:16]}  , {'Evidence', 'Position'}, {[], []}, 'towers', 'all', 'both');
        else
            disp('only works with Position or Evidence currently');
        end
        [~, topROI] = sort(out_Data.skaggsMetric.skaggs_real,'descend');
        dataUsed = dataUsed_full(:,topROI(1:round(length(topROI).*.1)));
    end
       
    if strcmp(regType,'GP')
        mdl = fitrgp(dataUsed(CV(i).trainLocations,:), trainResp1,'Standardize',false);%,'OptimizeHyperparameters','all');
        
        testResp_pred{i} = predict(mdl, dataUsed(CV(i).testLocations,:));
        trainResp_pred{i} = predict(mdl, dataUsed(CV(i).trainLocations,:));
        if exist('trainResp2')
            mdl2 = fitrgp(dataUsed(CV(i).trainLocations,:), trainResp2,'Standardize',false);%,'OptimizeHyperparameters','all');
            testResp_pred2{i} = predict(mdl2, dataUsed(CV(i).testLocations,:));
            trainResp_pred2{i} = predict(mdl2, dataUsed(CV(i).trainLocations,:));
        end
        
    elseif strcmp(regType,'SVM')
        mdl = fitrsvm(dataUsed(CV(i).trainLocations,:), trainResp1,'Standardize',true,'KernelFunction','gaussian');
        %mdl = fitrsvm(manifold3d(CV(i).trainLocations,:), trainResp1,'OptimizeHyperparameters','all');
        
        testResp_pred{i} = predict(mdl, dataUsed(CV(i).testLocations,:));
        trainResp_pred{i} = predict(mdl, dataUsed(CV(i).trainLocations,:));
        
    elseif strcmp(regType,'LLE')
        map1 = LLEMap([10:10:40],[0.00001 .001 1]);
        map1.fit(dataUsed(CV(i).trainLocations,:),trainResp1);
        testResp_pred{i} = map1.transform(dataUsed(CV(i).testLocations,:));
        trainResp_pred{i} = map1.transform(dataUsed(CV(i).trainLocations,:));
        
    elseif strcmp(regType,'GLM')
        mdl = fitglm(dataUsed(CV(i).trainLocations,:), trainResp1);%,'OptimizeHyperparameters','all');
        
        testResp_pred{i} = predict(mdl, dataUsed(CV(i).testLocations,:));
        trainResp_pred{i} = predict(mdl, dataUsed(CV(i).trainLocations,:));
    end
    
    testResp{i} = testResp1;
    
    if ~exist('trainResp2')
        trainResp{i} = trainResp1;
        
        nanind = find(isnan(testResp{i}));
        predLocations = CV(i).testLocations;
        if ~isempty(nanind)
            testResp{i}(nanind)=[];
            testResp_pred{i}(nanind)=[];
            predLocations(nanind) = [];
        end
        corr2(i) = corr(testResp{i},testResp_pred{i});
        dist2(i) = sqrt(mean((testResp{i}-testResp_pred{i}).^2));
    else
        trainResp{i} = trainResp1 - trainResp2;
        corr2(i) = corr(testResp{i},testResp_pred{i}-testResp_pred2{i});
        dist2(i) = sqrt(mean((testResp{i}-(testResp_pred{i}-testResp_pred2{i})).^2));
    end
    disp(['Cross validation fold ' num2str(i) ' of ' num2str(numFolds)]);
    if manifoldToggle~=1 && isempty(roiList)
        topROIlist{i} = topROI;
    end
    
    %% For binary variables
    if strcmp(varType, 'PriorChoice') || strcmp(varType, 'PriorCorrect') || strcmp(varType, 'Choice') || strcmp(varType, 'ChoiceCorrect')
        
        if strcmp(varType, 'PriorChoice')
            testRespBin = CV(i).testData.PriorChoice;
        elseif strcmp(varType, 'PriorCorrect')
            testRespBin = CV(i).testData.PriorCorrect;
        elseif strcmp(varType, 'Choice')
            testRespBin = CV(i).testData.Choice;
        elseif strcmp(varType, 'ChoiceCorrect')
            testRespBin = CV(i).testData.ChoiceCorrect;
        end
        
        
        testTrial1  = CV(i).testData.Trial;
        testResp_pred1 = testResp_pred{i};
        trialn = unique(testTrial1);
        for j=1:length(trialn)
           
            curTrial = trialn(j);
            testResp_pred_cur = testResp_pred1(testTrial1==curTrial);
            testResp_pos_cur  = posResp1(testTrial1==curTrial);
            testResp_cur = testRespBin(testTrial1==curTrial);
            resp1_pred = mean(testResp_pred_cur);
            resp1 = mean(testResp_cur);
            
            respAll(j,1) = resp1_pred;
            respAll(j,2) = resp1;
            respAll(j,3) = abs(resp1-resp1_pred);
            respAll(j,4) = abs(resp1-resp1_pred)<.5;
            
            respPos{i}{j} = testResp_pos_cur;
            respPred{i}{j} = testResp_pred_cur;
        end
        corPer(i) = sum(respAll(:,4))/size(respAll(:,4),1);
        corPer_raw{i} = respAll;
        clear respAll;
        posVals = vertcat(respPos{i}{:});
        respVals = vertcat(respPred{i}{:});
        
        [~,~,bin1] = histcounts(posVals,[-5:10:305]);
        for k=1:max(bin1)
           curBin = respVals(bin1==k); 
           binValmean(i,k) = mean(curBin);
        end
    end
    
    trueResponse(CV(i).testLocations)=testResp1;
    predResponse(predLocations)=testResp_pred{i};
    
end

% Put the data back into the same length as pre-manifold got-rid-of data
Datarange = outMind.Datarange;
reconData = zeros(size(Datarange));
realData = zeros(size(Datarange));
realData(Datarange) = trueResponse;
reconData(Datarange) = predResponse;

%% Save the outputs

outputNonlinearDecoding.corr2 = corr2;
outputNonlinearDecoding.meancorr = mean(corr2);
outputNonlinearDecoding.dist2 = dist2;
outputNonlinearDecoding.meandist = mean(dist2);
outputNonlinearDecoding.testResp = testResp;
outputNonlinearDecoding.testResp_pred = testResp_pred;
outputNonlinearDecoding.trainResp = trainResp;
outputNonlinearDecoding.trainResp_pred = trainResp_pred;
outputNonlinearDecoding.trueResponse = trueResponse;
outputNonlinearDecoding.predResponse = predResponse;
outputNonlinearDecoding.realData = realData;
outputNonlinearDecoding.reconData = reconData;
outputNonlinearDecoding.behavioralVariables = behavioralVariables;

if manifoldToggle~=1 && isempty(roiList)
    outputNonlinearDecoding.topROIlist = topROIlist;
end

if strcmp(varType, 'PriorChoice') || strcmp(varType, 'PriorCorrect') || strcmp(varType, 'Choice') || strcmp(varType, 'ChoiceCorrect')
   outputNonlinearDecoding.corPer = corPer;
   outputNonlinearDecoding.corPer_raw = corPer_raw;
   outputNonlinearDecoding.meancorPer = mean(corPer);
   outputNonlinearDecoding.respPos = respPos;
   outputNonlinearDecoding.respPred = respPred;
   outputNonlinearDecoding.binValmean = binValmean;
   outputNonlinearDecoding.binValmean_CV = mean(binValmean);
end