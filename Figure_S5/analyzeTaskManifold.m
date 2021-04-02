%close all;
%clear all;

%load('/Users/ms81/TT_embedding_manifold.mat')

%% Estimate dimensionality
data = result.forestdat.rwd.Dg;
xInputs = logspace(log10(0.0002),log10(0.4), 1000);
Cs = zeros(length(xInputs),length(data));
for p_idx = 1:length(data)
    distances = data(p_idx,:);
    distances(p_idx) = inf; % ignore distance to itself
    for idx = 1:length(xInputs)
       Cs(idx,p_idx) = sum( distances<xInputs(idx) );
    end
end

figure(1)
loglog(xInputs, mean(Cs,2),'.')
hold on;
xx = log(xInputs);
yy = log(mean(Cs,2)');
range = yy > -1 & yy < 5.5;
c = polyfit( xx(range), yy(range), 1);
loglog(xInputs, exp(c(2))* xInputs .^c(1), 'k--')
title(num2str(c(1)))
xlim(exp([min(xx(isfinite(yy))), max(xx(isfinite(yy)))]))
ylim(exp([min(yy(isfinite(yy))), max(yy(isfinite(yy)))]))
ylim([1e-1, 1e3])
xlim([0.02 , 0.5])
set(gcf,'color','w');
xlabel("distance [a.u.]")
ylabel("# Neighbors")


%% How well does MIND PCA do, compared to PCA?

testd = ismember(data_trial, 202:210);
test_flag = (data_phase(testd)<450) & (data_positionY(testd)>0);
test_data = all_test_data(test_flag, :); 

[coeff,score,latent,tsquared,explained,mu] = pca(data_for_mind); 
% to vet variance explained for PCA

cc = [];

for dim = 1:length(dembed)
    disp(dim)
    % check MIND
    pca_coords = result.forestdat.pca.model.transform(test_data, 0.95);
    y_testdata = result.allembed(dim).f2m.map.transform(pca_coords);
    pca_reconstructed = result.allembed(dim).m2f.map.transform(y_testdata);
    reconstructeed_test_data = result.forestdat.pca.model.inverse_transform(pca_reconstructed);
    cc = [cc; corr(test_data(:), reconstructeed_test_data(:))];
    
%     figure(dim)
%     subplot(2,1,1)
%     imagesc(reshape(all_test_data(100,:),length(dx),length(dy)))
%     title("Full frame 100")
%     subplot(2,1,2)
%     imagesc(reshape(reconstructeed_test_data(100,:),length(dx),length(dy)))
%     title(['Reconstructed from d=', num2str(dim)])
end

figure(2)
bar([cc.^2, cumsum(explained(1:length(cc))/100)])
title("MIND vs. PCA - R^2")
ylabel("crossvalidated R^2")
legend(["MIND", "PCA"])
ylim([0,0.7])
xlabel("# of dimensions")
set(gcf,'color','w');


%% And the manifold movies
% flag is the relevant range for fit: (data_phase<450) & (data_positionY>0)
pca_coords = result.forestdat.pca.model.transform(all_data(flag, :), 0.95);
y = result.allembed(3).f2m.map.transform(pca_coords);
    
% make_movie(y, [], 20, data_positionX(flag), 't_positionX.mp4')
% make_movie(y, [], 20, data_positionY(flag), 't_positionY.mp4')
% make_movie(y, [], 20, data_phase(flag), 't_phase.mp4')
% make_movie(y, [], 20, data_luminance(flag), 't_luminance.mp4')
% make_movie(y, [], 20, data_rightT(flag) - data_leftT(flag), 't_deltaT.mp4')
% make_movie(y, [], 20, mod(data_viewangle(flag)+pi,2*pi)-pi, 't_viewangle.mp4')

plot_manifold_figure(y, data_evidence(flag), 45, 35, -9.85, 7.55, '/Users/ms81/', true)
plot_manifold_figure(y, data_luminance(flag), 45, 35, 0, 30, '/Users/ms81/', true)


