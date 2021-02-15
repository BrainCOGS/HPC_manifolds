function outputPlotHPC2HPC = plot_HPC2HPC(outputHPC2HPC)

% Data from the best crossvalidated GPR:
GPR_evidence = diag(outputHPC2HPC.M_evidence)';
GPR_position = diag(outputHPC2HPC.M_position)';

M_evidence = outputHPC2HPC.M_evidence - diag(GPR_evidence);
M_position = outputHPC2HPC.M_position - diag(GPR_position);


bar_evi = [];
bar_Y = [];
for ani = 1:7
    dataE = M_evidence(:,ani);
    bar_evi = [bar_evi, max(dataE(dataE~=0))];
    dataY = M_position(:,ani);
    bar_Y = [bar_Y, max(dataY(dataY~=0))];
end

%[h,pE] = ttest(bar_evi, GPR_evidence);
%[h,pY] = ttest(bar_Y, GPR_position);

pE = signrank(bar_evi, GPR_evidence);
pY = signrank(bar_Y, GPR_position);

figure;
subplot(1,3,1)
bar([mean(GPR_evidence), mean(bar_evi)])
hold on;
plot([GPR_evidence; bar_evi], 'ko-')
errorbar([1, 2], [mean(GPR_evidence), mean(bar_evi)], [std(GPR_evidence)/sqrt(length(GPR_evidence)), std(bar_evi)/sqrt(length(bar_evi))], 'bo');
title("evidence: p=" + num2str(pE))
xlabel("GPR, Hyper")
ylabel("Reconstruction Score (r) for Evidence")

subplot(1,3,2)
bar([mean(GPR_position), mean(bar_Y)])
hold on;
plot([GPR_position; bar_Y], 'ko-')
errorbar([1, 2], [mean(GPR_position), mean(bar_Y)], [std(GPR_position)/sqrt(length(GPR_position)), std(bar_Y)/sqrt(length(bar_Y))], 'bo');
title("position: p=" + num2str(pY))
xlabel("GPR, Hyper")
ylabel("Reconstruction Score (r) for Position")

subplot(1,3,3)
dY = bar_Y.^2 ./ max([bar_Y.^2 ; GPR_position.^2]);
dE = bar_evi.^2 ./ max([bar_evi.^2 ; GPR_evidence.^2]);
%[h,pYvE] = ttest(dY, dE);
pYvE = signrank(dY, dE);
bar([mean(dY), mean(dE)])
hold on;
plot([dY; dE], 'ko-')
errorbar([1, 2], [mean(dY), mean(dE)], [std(dY)/sqrt(length(dY)), std(dE)/sqrt(length(dE))], 'bo');
title("p=" + num2str(pYvE))
xlabel("Pos-shared, Evi-shared")
ylabel("Contribution of shared geometry to Reconstruction score")
set(gcf,'color','w');

outputPlotHPC2HPC.pE = pE;
outputPlotHPC2HPC.pY = pY;
outputPlotHPC2HPC.pYvE = pYvE;
outputPlotHPC2HPC.dY_mean = mean(dY);
outputPlotHPC2HPC.dY_sem  = nieh_sem(dY);
outputPlotHPC2HPC.dE_mean = mean(dE);
outputPlotHPC2HPC.dE_sem  = nieh_sem(dE);

