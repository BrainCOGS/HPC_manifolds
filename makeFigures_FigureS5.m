
%% Run the code

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Run this code to generate the outMind variable (task manifold)
% dataPath = 'C:\Neuroscience\imaging\FINAL\taskmanifold_Data\TT_embedding_rawdata.mat';
% saveName = 'TT_embedding_manifold.mat';
% fitManifold2Video(fnameStruct(7).fname, '', dataPath, saveName, false);
% Instead, to keep from rerunning manifold fitting, load data:
load("C:\Neuroscience\imaging\FINAL\taskmanifold_Data\TT_embedding_manifold.mat")

% Run this code to get the luminance values 
load('C:\Neuroscience\imaging\FINAL\taskmanifold_Data\luminace2HPC.mat')
% Or run this code:
% fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');
% videopath = '\\192.168.0.233\Neuroscience\schottdorf\E65_trials\';
% vfilename = 'hnieh_E65-2018-02-02-B';
% outputLuminance2HPC = luminance2HPC(fnameStruct, videopath, vfilename);

luminance = outputLuminance2HPC.luminance;
trial_luminance = outputLuminance2HPC.trial_luminance;

%% Make the luminance plot

figure;
flag = (data_phase<450) & (data_positionY>0);
scatter3(data_positionX(flag), data_positionY(flag), data_luminance(flag), 5, 'filled','MarkerFaceAlpha',.5)
sourceData_S5a = [double(data_positionX(flag)) double(data_positionY(flag)) data_luminance(flag)];
view(58.5, 30);
set(gcf, 'Position', [500, 340, 660, 420]);
set(gcf,'renderer','painters');


%% Make the plots for the task manifold

pca_coords = result.forestdat.pca.model.transform(data_for_mind, 0.95);
y = result.allembed(3).f2m.map.transform(pca_coords);

smooth_DL = mind_smoothDimensions(y,data_luminance(flag),20);
smooth_DE = mind_smoothDimensions(y,data_evidence(flag),20);

load(fnameStruct(7).fname_mani);
dimEmbed = 3;
sample = outMind.dataDFF;
pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
manifold3d = outMind.dat.allembed(outMind.config_input.mindparameters.embed.d==dimEmbed).f2m.map.transform(pca_coords);

lum = trial_luminance(outMind.Datarange);
smooth_L = mind_smoothDimensions(manifold3d,lum,20);

figure;
subplot(1,3,1)
scatter3(y(:,1), y(:,2), y(:,3), 20, smooth_DL,'filled','MarkerFaceAlpha',.5)
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')
title('Luminance (smoothed)')
xlim([-0.14,0.14])
ylim([-0.14,0.14])
zlim([-0.14,0.14])
axis square
colorbar
grid on;
view(118,29)
caxis([0 30])

subplot(1,3,2)
scatter3(y(:,1), y(:,2), y(:,3), 20, smooth_DE,'filled','MarkerFaceAlpha',.5)
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')
title('Evidence (smoothed)')
xlim([-0.14,0.14])
ylim([-0.14,0.14])
zlim([-0.14,0.14])
axis square
colorbar
grid on;
view(118,29)
% Taken from the color axis from Figure 3
caxis([-9.85 7.55])

sourceData_S5b = [y(:,1) y(:,2) y(:,3) smooth_DL' smooth_DE'];

% Make plots for the neural manifold
subplot(1,3,3)
scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3), 20, smooth_L,'filled','MarkerFaceAlpha',.5);
xlabel('Dim 1')
ylabel('Dim 2')
zlabel('Dim 3')
title('Luminance (smoothed)')
xlim([-0.07,0.07])
ylim([-0.07,0.07])
zlim([-0.07,0.07])
axis square
colorbar
grid on;
view(-63,36)
caxis([0 30])

sourceData_S5c = [manifold3d(:,1) manifold3d(:,2) manifold3d(:,3) smooth_L'];
% The evidence plot is identical to the one in figure 3

set(gcf, 'Position', [100, 340, 1300, 420])
set(gcf,'renderer','painters');
