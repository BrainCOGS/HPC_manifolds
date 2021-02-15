function outputReconstructionMultipleROIs = collect_allTrials_reconstructionMultipleROIs_Function(fnameStruct, fileLoc, aniN, roiList)

% fileLoc is the folder where the files are, i.e. D:\CrossValidation_holdOneCell\
% aniN is the name of the animal, i.e. 'E43' or a number 1-7
% roiList is the held-out rois to use, if blank, will look through folder
% and do all of them.

%% Set up variables
argins.fnameStruct = fnameStruct;
argins.fileLoc     = fileLoc;
argins.aniN        = aniN;
argins.roiList     = roiList;

outputReconstructionMultipleROIs.argins = argins;

config.inputRNG = 1;             % random seed
rng(config.inputRNG);
nleafs       = 500;              % list of leafs for fitting manifolds of different complexity
lmf          = 1;                % list of landmark_fractions
dimList      = [2:7];
% Only do the 5 dim manifold if we're plotting the reconstruction figure panel
if length(roiList)==40
    dimList = 5;
end
config.nleafs = nleafs;
config.lmf = lmf;
config.dimList = dimList;
outputReconstructionMultipleROIs.config = config;

%% if aniN is a number, make it the right animal name
if isnumeric(aniN)
    if aniN==1
        aniN = 'E22';
    elseif aniN==2
        aniN = 'E39';
    elseif aniN==3
        aniN = 'E43';
    elseif aniN==4
        aniN = 'E44';
    elseif aniN==5
        aniN = 'E47';
    elseif aniN==6
        aniN = 'E48';
    elseif aniN==7
        aniN = 'E65';
    end
end

%% Update file location based on animal

fileLoc = [fileLoc aniN '/allTrials/']

%% If roiList is empty, get all roiout folder from the directory

if isempty(roiList)
    listing = dir(fileLoc);
    names1 = [listing.name];
    roiList = strsplit(names1, 'roiout_');
    roiList = roiList(2:end);
    roiList = cellfun(@str2num, roiList);
end

%% Start analysis

if strcmp(aniN, 'E22')
    fname = fnameStruct(1).fname;
elseif strcmp(aniN, 'E39')
    fname = fnameStruct(2).fname;
elseif strcmp(aniN, 'E43')
    fname = fnameStruct(3).fname;
elseif strcmp(aniN, 'E44')
    fname = fnameStruct(4).fname;
elseif strcmp(aniN, 'E47')
    fname = fnameStruct(5).fname;
elseif strcmp(aniN, 'E48')
    fname = fnameStruct(6).fname;
elseif strcmp(aniN, 'E65')
    fname = fnameStruct(7).fname;
else
    error("no .modeling.mat file found")
end

nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [0 0], fname,'none','towers', 1, 1);
behavioralVariablesFull = nic_output.behavioralVariables;
ROIactivitiesFull = nic_output.ROIactivities;

for iroi = 1:length(roiList)
    for idim = 1:length(dimList);
        curDim = dimList(idim);
        leaveoutROI = roiList(iroi);
        folder_t = [fileLoc, 'roiout_', num2str(leaveoutROI), '/'];
        
        f = [folder_t, 'minLeaves_' num2str(nleafs), '_lmf_', num2str(lmf), '_crossval.mat'];
        f2 = [folder_t, 'minLeaves_' num2str(nleafs), '_lmf_', num2str(lmf), '_manifold.mat'];
        disp(f);
        
        if isfile(f)
            load(f);
            
            outputNonlinearDecoding_roiout = mind_nonlinearDecoding_dimX_All(fname, f2,5,'GP',leaveoutROI,'towers',0,1,[], curDim);
            reconstructedAll(:,iroi, idim) = outputNonlinearDecoding_roiout.reconData;
            corrAll(iroi, idim) = outputNonlinearDecoding_roiout.meancorr;
        end
        disp(['Finished dim ' num2str(curDim) ' in ' num2str(iroi) ' roi of ' num2str(length(roiList))]);
    end
end

outputReconstructionMultipleROIs.reconstructedAll = reconstructedAll;
outputReconstructionMultipleROIs.corrAll = corrAll;


%% Collect and plot if it's E43 and it's the full 40 rois

if strcmp(aniN, 'E43')==1 && length(roiList)==40
    
    % set to 1 less than desired dim
    % reconstructedAll = squeeze(reconstructedAll(:,:,4));
    
    trialsList = [145 146 147 148 149 151 152 153 154 157 158 159 160];
    nic_output = extractVariables('all', 2, trialsList, 2, 0, 0, 5, [0 0], fname,'none','towers',1,1);
    ROIactivities = nic_output.ROIactivities;
    
    nic_output2 = extractVariables('all', 2, trialsList, 2, 0, 0, 5, [11 4], fname,'none','towers',1,1);
    ROIactivities_11_4 = nic_output2.ROIactivities;
    
    ROIactivities_11_4_Datarange    = sum(ROIactivities_11_4,2)>0;
    ROIactivities_11_4 = ROIactivities_11_4(ROIactivities_11_4_Datarange,:);
    ROIactivities = ROIactivities(ROIactivities_11_4_Datarange,:);
    
    % Get the reconstructedAll data that's only in the desired trials
    trialsKeep = ismember(behavioralVariablesFull.Trial,trialsList);
    justTrialData = reconstructedAll(trialsKeep,:);
    justTrialData = justTrialData(ROIactivities_11_4_Datarange,:);
    reconstructedAll_11_4 = reconstructedAll - mean(reconstructedAll);
    reconstructedAll_11_4 = mind_preprocess(reconstructedAll_11_4,11,4);
    justTrialData_11_4 = reconstructedAll_11_4(trialsKeep,:);
    justTrialData_11_4 = justTrialData_11_4(ROIactivities_11_4_Datarange,:);
    
    mean(corrAll)
    std(corrAll)
    
    figure;
    ax1 = subplot(2,2,1);
    imagesc(mat2gray(ROIactivities(400:600,1:40)'));
    xticks([0 5*(1/.0694) 10*(1/.0694)])
    xticklabels({'0' '5' '10'});
    set(gca, 'box', 'off')
    ylabel('Neuron #');
    title('raw')
    axis square
    
    ax2 = subplot(2,2,2);
    imagesc(mat2gray(ROIactivities_11_4(400:600,1:40)'));
    xticks([0 5*(1/.0694) 10*(1/.0694)])
    xticklabels({'0' '5' '10'});
    set(gca, 'box', 'off')
    ylabel('Neuron #');
    title('raw (smoothed and thresholded)')
    axis square
    
    ax3 = subplot(2,2,3);
    imagesc(mat2gray(justTrialData(400:600,1:40)'));
    xticks([0 5*(1/.0694) 10*(1/.0694)])
    xticklabels({'0' '5' '10'});
    set(gca, 'box', 'off')
    xlabel('Time (s)');
    title('reconstructed')
    axis square
    
    ax4 = subplot(2,2,4);
    imagesc(mat2gray(justTrialData_11_4(400:600,1:40)'));
    xticks([0 5*(1/.0694) 10*(1/.0694)])
    xticklabels({'0' '5' '10'});
    set(gca, 'box', 'off')
    title('reconstructed (smoothed and thresholded)')
    axis square
    
    suptitle('Reconstruct Held-Out ROI');
    
end
