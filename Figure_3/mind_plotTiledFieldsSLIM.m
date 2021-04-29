function outputTiledFields = mind_plotTiledFieldsSLIM(fname, fname_mani)

% Load colormap
load('fieldwidth_colormap2.mat')

% Set up some config variables
numDim=3;
config.numDim = numDim;
stdThreshold = 3;
config.stdThreshold = stdThreshold;
load(fname_mani);

outputTiledFields.config = config;

%% Get manifold Data
sample = outMind.dataDFF;
pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
manifold3d = outMind.dat.allembed(outMind.config_input.mindparameters.embed.d==numDim).f2m.map.transform(pca_coords);
Datarange = outMind.Datarange;

%% Fit Firing Fields
outputFitFiringFields = mind_fitFiringFieldsNEW_dimX_manuel(fname, fname_mani,'all',[5 2],1,5,0,5,5,0,3);
ROIscores = outputFitFiringFields.meanCorr;
[sort1, sort2] = sort(ROIscores, 'descend');
keepROIs = sort2(1:25);
activeROIs = outputFitFiringFields.activeROIsWidth;
outputTiledFields.keepROIs = keepROIs;

nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [5 2], fname,'none','towers',1,1);
ROIactivities = nic_output.ROIactivities;
ROIactivities = ROIactivities(:,activeROIs);
ROIactivities = ROIactivities(:,keepROIs);
ROIactivities = ROIactivities(Datarange,:);


%% Using just the final ROIactivities, find points above a number of std

Data_mean = nanmean(ROIactivities);
Data_std  = nanstd(ROIactivities);
Data_thres = Data_mean + (Data_std.*stdThreshold);
Data_thres = repmat(Data_thres,size(ROIactivities,1),1);

tempThres = ROIactivities>Data_thres;
ROIactivities_thres = ROIactivities.*tempThres;


%% Plot a single example and 5 on top of each other

plotROI = [3 14 25 18 9];

colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330];

figure;
subplot(1,2,1)
curROI = plotROI(2);
hold on;
ROInonzero = find(ROIactivities_thres(:,curROI));
ROIzero    = boolean(ones(length(ROIactivities_thres),1));
ROIzero(ROInonzero) = 0;
h          = scatter3(manifold3d(ROInonzero,1),manifold3d(ROInonzero,2),manifold3d(ROInonzero,3), 20, ROIactivities_thres(ROInonzero,curROI),'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8]);
sourceData_3d_cell = [manifold3d(ROInonzero,1) manifold3d(ROInonzero,2) manifold3d(ROInonzero,3)      ROIactivities_thres(ROInonzero,curROI)];
h2         = scatter3(manifold3d(ROIzero,1),manifold3d(ROIzero,2),manifold3d(ROIzero,3), 20, 'filled', 'MarkerFaceColor', [.8 .8 .8], 'MarkerFaceAlpha', .05);
sourceData_3d_grey = [manifold3d(ROIzero,1) manifold3d(ROIzero,2) manifold3d(ROIzero,3)];
colormap(cmap_firingfield);
caxis([0 2]);
xlim([-0.07,0.07])
ylim([-0.07,0.07])
zlim([-0.07,0.07])
%view(-51,-35)
view(-63,36)
%axis tight
axis square
grid on;
xlabel('Dim 1');
ylabel('Dim 2');
zlabel('Dim 3');

subplot(1,2,2)
ROIzeroALL = ones(size(ROIactivities,1),1);
ROIzeroALL = boolean(ROIzeroALL);
for i=1:length(plotROI)
    curROI = plotROI(i);
    hold on;
    ROInonzero = find(ROIactivities_thres(:,curROI));
    ROIzeroALL(ROInonzero) = 0;
    h             = scatter3(manifold3d(ROInonzero,1),manifold3d(ROInonzero,2),manifold3d(ROInonzero,3), 20, colors(i,:),'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8]);
    sourceData_3e_cell{i} = [manifold3d(ROInonzero,1) manifold3d(ROInonzero,2) manifold3d(ROInonzero,3)];
end
h2         = scatter3(manifold3d(ROIzeroALL,1),manifold3d(ROIzeroALL,2),manifold3d(ROIzeroALL,3), 20, 'filled', 'MarkerFaceColor', [.8 .8 .8], 'MarkerFaceAlpha', .05);
sourceData_3e_grey = [manifold3d(ROIzeroALL,1) manifold3d(ROIzeroALL,2) manifold3d(ROIzeroALL,3)];
colormap(cmap_firingfield);
caxis([0 2]);
xlim([-0.07,0.07])
ylim([-0.07,0.07])
zlim([-0.07,0.07])
view(-63,36)
axis square
grid on;
xlabel('Dim 1');
ylabel('Dim 2');
zlabel('Dim 3');
set(gcf,'renderer','painters');


%% Save the variables

outputTiledFields.manifold3d = manifold3d;
outputTiledFields.ROIactivities_thres = ROIactivities_thres;
outputTiledFields.sourceData_3d_cell = sourceData_3d_cell;
outputTiledFields.sourceData_3d_grey = sourceData_3d_grey;
outputTiledFields.sourceData_3e_cell = sourceData_3e_cell;
outputTiledFields.sourceData_3e_grey = sourceData_3e_grey;