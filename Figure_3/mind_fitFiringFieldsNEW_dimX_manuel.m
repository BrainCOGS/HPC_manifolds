function outputFitFiringFields = mind_fitFiringFieldsNEW_dimX_manuel(fname, fname_mani, whichROIs, preprocessParam, toggleActive, eventThreshold, toggleSpock, numFolds, numFoldsEval, toggleShuffle, dimManifold)
% toggleActive == 1 means take only cells with transients in more than 10%
% trials
% eventThreshold is used for detecting active cells, either 11 or 5
% toggleSpock changes code for spock or laptop
% toggleShuffle, shuffles the data if equal to 1
rng(1);

if toggleSpock==1
    pc = parcluster('local')
    
    pc.JobStorageLocation = strcat('/tmp/',getenv('USER'),'-',getenv('SLURM_JOB_ID'))
    
    parpool(pc,12)
    
    p=gcp;
    
    poolSize = p.NumWorkers
    
end

tic
argins.fname = fname;
argins.fname_mani = fname_mani;
argins.whichROIs = whichROIs;
argins.preprocessParam = preprocessParam;
argins.toggleActive    = toggleActive;
argins.eventThreshold  = eventThreshold;
argins.numFolds        = numFolds;
argins.numFoldsEval    = numFoldsEval;
argins.toggleShuffle   = toggleShuffle;
argins.dimManifold     = dimManifold;

load(fname_mani)
load(fname)

numDim=dimManifold;
sample = outMind.dataDFF;
pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
manifold3d = outMind.dat.allembed(outMind.config_input.mindparameters.embed.d==numDim).f2m.map.transform(pca_coords);
Datarange = outMind.Datarange;

nic_output = extractVariables(whichROIs, 2, 'keepTrials', 2, 0, 0, 5, preprocessParam, fname,'none','towers',1,1);
behavioralVariables = nic_output.behavioralVariables;
behavioralVariables = behavioralVariables(Datarange,:);

trialn   = behavioralVariables.Trial;
trials = unique(trialn);

ROIactivities = nic_output.ROIactivities;
dataDFF = ROIactivities(Datarange,:);

if toggleActive==1
    % Using [5 0] as preprocessParam because that's what we defined an
    % event as for doublets
    nic_output2 = extractVariables(whichROIs, 2, 'keepTrials', 1, 0, 0, eventThreshold, [5 0], fname,'none','towers',1,1);
    ROIactivities2 = nic_output2.ROIactivities;
    numtransients = sum(ROIactivities2);
    activeROIsWidth = numtransients>(length(unique(trialn))*.1);
    dataDFF = dataDFF(:,activeROIsWidth);
end

if toggleShuffle == 1
    jeff_constant = 10;
    dataDFF = jeff_randomizer(dataDFF, jeff_constant);
end


% outputDoublets_E65_new5 = findDoublets_20190723(5, 0, 0, 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat',  [0 0],1);
% cell1 = dataDFF(:,4);
% cell2 = dataDFF(:,12);
% doubTrials = outputDoublets_E65_new5.saveAll_basics.aTrials(outputDoublets_E65_new5.saveAll_doublets_Sig(16).trials_appear);

% trialNums = behavioralVariables.Trial;
% subTrials = ismember(trialNums, doubTrials);
%
% behavioralVariables = behavioralVariables(subTrials,:);
% dataDFF=dataDFF(subTrials,:);
% manifold3d = manifold3d(subTrials,:);


%%
%numFolds = 5;
CV = generateCrossValSet_v2(behavioralVariables, numFolds);
opt_k3 = zeros(size(dataDFF,2),numFoldsEval);
opt_lambda3 = zeros(size(dataDFF,2),numFoldsEval);
%predDFF = zeros(size(dataDFF));
%corrPred = zeros(size(dataDFF,2),1);
%if strcmp(whichROIs,'all')
whichROIs = size(dataDFF,2);
%end


%parfor j=1:whichROIs%1:size(dataDFF,2)
for j=1:whichROIs%1:size(dataDFF,2)
    
    cell1 = dataDFF(:,j);
    
%     CV = allCV(j);
%     CV = CV.CV;
    
    for i=1:numFoldsEval
        
        trainData3 = manifold3d(CV(i).trainLocations,:);
        testData3 = manifold3d(CV(i).testLocations,:);
        
        trainLabels_1 = cell1(CV(i).trainLocations);
        testLabels_1  = cell1(CV(i).testLocations);
        
        if sum(trainLabels_1>0)>2 && sum(testLabels_1>0)>2
            
            win = trainLabels_1>0;
            %fo = min(trainLabels_1(win));
            
            win2 = testLabels_1>0;
            %fo2 = min(testLabels_1(win2));

            trainLabels_1 = trainLabels_1(win);
            testLabels_1  = testLabels_1(win2);
            
            trainData3 = trainData3(win,:);
             
            testData3 = testData3(win2,:);
            
            %map1 = LLEMap([12:5:50],[10.^(-8:1:0)]);
            if length(trainData3)<50
                disp('in if statement where not enough trainData');
                map3 = LLEMap([10:10:(length(trainData3)-5)],[0.00001 1]);
            else
                map3 = LLEMap([10:10:40],[0.00001 1]);
            end
            %map1 = LLEMap([10 40],[0.00001 1]);
            map3.fit(trainData3,trainLabels_1);
            opt_k3(j,i) = map3.k;
            opt_lambda3(j,i) = map3.lambda;
            testResp_pred_3  = map3.transform(testData3);
            
            % corr_1(i) = corr(testResp_pred_1,testLabels_1);
            % corr_2(i) = corr(testResp_pred_2,testLabels_2);
            
            %predDFF(CV(i).testLocations,j) = testResp_pred_1;
            % cell2_pred(CV(i).testLocations) = testResp_pred_2;
            corrCV_3d(j,i) = corr(testResp_pred_3,testLabels_1);
             
            dist_3d(j,i)  = mean(abs((mat2gray(testResp_pred_3)-mat2gray(testLabels_1))));
             
            pred3d{j,i} = testResp_pred_3;
            
        else
            
            corrCV_3d(j,i) = NaN;
            
            dist_3d(j,i)  = NaN;
             
            pred3d{j,i} = NaN;
            
        end
    end
    %corrPred(j) = corr(predDFF(:,j),dataDFF(:,j));
    j
end

%
% figure;
% subplot(2,2,1);
% scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3),5, cell1);
% subplot(2,2,2);
% scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3),5, cell2);
% subplot(2,2,3);
% scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3),5, cell1_pred);
% subplot(2,2,4);
% scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3),5, cell2_pred);


%outputDoubletTrajectory.predDFF = predDFF;
%outputDoubletTrajectory.corrPred = corrPred;
outputFitFiringFields.opt_k    = opt_k3
outputFitFiringFields.opt_lambda = opt_lambda3;
outputFitFiringFields.corrCV = corrCV_3d;
outputFitFiringFields.meanCorr = nanmean(corrCV_3d,2);
outputFitFiringFields.meanDist = nanmean(dist_3d,2);
outputFitFiringFields.allmeanCorr = mean(nanmean(corrCV_3d,2));
outputFitFiringFields.allmeanDist = mean(nanmean(dist_3d,2));
outputFitFiringFields.pred  = pred3d;
outputFitFiringFields.argins     = argins;
outputFitFiringFields.activeROIsWidth = activeROIsWidth;
outputFitFiringFields.manifold = manifold3d;

elapsedTime = toc
outputFitFiringFields.elapsedTime = elapsedTime;
