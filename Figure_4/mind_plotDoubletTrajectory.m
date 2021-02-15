function [outputTrajectory] = mind_plotDoubletTrajectory(fname_mani, fname, togglePlots, doubletSD, preprocessparam, numDim, outputDoublets, doubletNum, selectedTrials, sigToggle)

% if outputDoublets is empty, it'll run the find doublets function
% if sigToggle==1, then use the significantly predictive doublets
% If togglePlots==1, make the plots, otherwise, skip
% preprocessparam is ONLY for the doublet finding
% 
% Sample Call: 
% [outputTrajectory, config] = mind_plotTrajectory('E47_activeROIs_removeITI_onlyGoodMain_10628dp_219ROI_100trees_5314land_dim2-6_11pre_noThreshold.mat', 'C:\Neuroscience\imaging\E47\20170927\E47_20170927_70per_userSetSD5minDur0.modeling.mat', togglePlots);

%% Save inputs and set up config

argins.fname_mani = fname_mani;
argins.fname = fname;
argins.togglePlots = togglePlots;
argins.doubletSD = doubletSD;
argins.preprocessparam = preprocessparam;
argins.numDim = numDim;
argins.outputDoublets = outputDoublets;
argins.doubletNum = doubletNum;
argins.selectedTrials = selectedTrials;
argins.sigToggle = sigToggle;

outputTrajectory.argins = argins;

config.rngInput = 42;
rng(config.rngInput)

outputTrajectory.config = config;


%% Get the doublets and manifold data

load(fname);
load(fname_mani);

sample = outMind.dataDFF;
pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.dat.mindparameters.pca.n);
manifold3d = outMind.dat.allembed(outMind.dat.mindparameters.embed.d==numDim).f2m.map.transform(pca_coords);

padSpots = find(outMind.Datarange==0);
manifold3d_pad = zeros(length(outMind.Datarange),numDim);
manifold3d_pad(outMind.Datarange,:) = manifold3d;
for i=1:length(padSpots)
   manifold3d_pad(padSpots(i),:) = manifold3d_pad(padSpots(i)-1,:); 
end
manifold3d = manifold3d_pad; 

nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [11 4], fname,'none','towers',1,1);
behavioralVariables = nic_output.behavioralVariables;

% Resample the positions so it's the same as the manifold dimensions
behavioralVariables = behavioralVariables(outMind.Datarange,:);
Positions = behavioralVariables.Position;
Positions_pad = zeros(length(outMind.Datarange),1);
Positions_pad(outMind.Datarange,:) = Positions;
for i=1:length(padSpots)
   Positions_pad(padSpots(i),:) = Positions_pad(padSpots(i)-1,:); 
end
Positions = Positions_pad;
behavioralVariables = nic_output.behavioralVariables;
Times = behavioralVariables.Time;

% Normalizes max distance to 1
Positions = mat2gray(Positions);
manifold3d(:,1) = mat2gray(manifold3d(:,1));
manifold3d(:,2) = mat2gray(manifold3d(:,2));
manifold3d(:,3) = mat2gray(manifold3d(:,3));
manifold3d = manifold3d.*(sqrt(1/3));

% Get the doublets
if isempty(outputDoublets)
    outputDoublets = findDoublets_20190509(doubletSD, 0, 0, fname,preprocessparam,1);
end

% Collect some variables
saveAll_doublets        = outputDoublets.saveAll_doublets;
saveAll_doublets_Sig    = outputDoublets.saveAll_doublets_Sig;
saveAll_preprocess      = outputDoublets.saveAll_preprocess;

outputTrajectory.fn_modeling = fname;
outputTrajectory.saveAll_doublets = saveAll_doublets;
outputTrajectory.saveAll_doublets_Sig = saveAll_doublets_Sig;
outputTrajectory.saveAll_preprocess = saveAll_preprocess;

trialNums = behavioralVariables.Trial;
trialList = unique(trialNums);

for i=1:length(trialList)
    manifold_splitted_all{i} = manifold3d(trialNums==trialList(i),:);
    Position_splitted_all{i} = Positions(trialNums==trialList(i));
    time_splitted_all{i}     = Times(trialNums==trialList(i));
end

outputTrajectory.manifold_splitted_all = manifold_splitted_all;
outputTrajectory.Position_splitted_all = Position_splitted_all;
outputTrajectory.time_splitted_all     = time_splitted_all;



%% Set up some variables to be used by both plotting trajectories and positions

outputTrajectory.config.exampleDoublet = doubletNum;

if sigToggle==1
    cell1 = saveAll_doublets_Sig(doubletNum).first_cell;
    cell2 = saveAll_doublets_Sig(doubletNum).second_cell;
    trials_appear = saveAll_doublets_Sig(doubletNum).trials_appear;
else
    cell1 = saveAll_doublets(doubletNum).first_cell;
    cell2 = saveAll_doublets(doubletNum).second_cell;
    trials_appear = saveAll_doublets(doubletNum).trials_appear;
end
saveAll_preprocess_trials = saveAll_preprocess(trials_appear);
manifold_splitted = manifold_splitted_all(trials_appear);
position_splitted = Position_splitted_all(trials_appear);
time_splitted     = time_splitted_all(trials_appear);

colorchoices = rand(100,3);

for i=1:length(trials_appear)
    trialActivity{i} = saveAll_preprocess_trials(i).digitized_eventStart(:,[cell1 cell2]);
    lengths1(i) = size(manifold_splitted{i},1);
end


%% Plot the transients of a doublet on the manifold trajectory in 3D

if togglePlots==1
    figure
    for j=1:length(trials_appear)
        
        curManifold = manifold_splitted{j};
        curTrial = trialActivity{j};
        curTrial = logical(curTrial);
        hold on;
        scatter3(curManifold(curTrial(:,1),1), curManifold(curTrial(:,1),2), curManifold(curTrial(:,1),3),'b','filled','MarkerFaceAlpha',.5);
        scatter3(curManifold(curTrial(:,2),1), curManifold(curTrial(:,2),2), curManifold(curTrial(:,2),3),'r','filled','MarkerFaceAlpha',.5);
        
        if j==selectedTrials(1) || j==selectedTrials(2)
            plot3(manifold_splitted{j}(find(curTrial(:,1)):find(curTrial(:,2)),1), manifold_splitted{j}(find(curTrial(:,1)):find(curTrial(:,2)),2), manifold_splitted{j}(find(curTrial(:,1)):find(curTrial(:,2)),3))
        end
    end
end
axis square
grid on;
xlabel('dim 1')
ylabel('dim 2')
zlabel('dim 3')
view(-222.1424, 10.8865)
title(['doubletnum: ' num2str(doubletNum) ', blue cell 1, red cell2']);

