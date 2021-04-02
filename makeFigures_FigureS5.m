
%% Run the code

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Run this code to generate the outMind variable (task manifold)
% dataPath = 'C:\Neuroscience\imaging\FINAL\taskmanifold_Data\TT_embedding_rawdata.mat';
% saveName = 'TT_embedding_manifold.mat';
% fitManifold2Video(fnameStruct(7).fname, '', dataPath, saveName, false);
% Or load the data:
load('\\192.168.0.233\Neuroscience\schottdorf\fit_manifold2videodata\TT_embedding_manifold.mat')

% Run this code to get the luminance values 
luminance2HPC;

% Make the plots for the task manifold
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
scatter3(y(:,1), y(:,2), y(:,3), [], smooth_DL, '.')
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
scatter3(y(:,1), y(:,2), y(:,3), [], smooth_DE, '.')
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

set(gcf, 'Position', [100, 340, 1300, 420])

