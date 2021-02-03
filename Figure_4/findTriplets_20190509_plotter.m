function findTriplets_20190509_plotter(outputTriplets)

rng(1)

figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,3,1)
plot([outputTriplets.saveAll_triplets_Sig_Left.triplet_shuffle],'LineWidth',2);
hold on;
plot([outputTriplets.saveAll_triplets_Sig_Left.prediction],'LineWidth',2);
legend('Scamble','Real');
title('Significant Left-Preferring Triplets');
xlabel('Triplet');
ylabel('Fraction Trials Going Left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,2)
plot([outputTriplets.saveAll_triplets_Sig_Right.triplet_shuffle],'LineWidth',2);
hold on;
plot([outputTriplets.saveAll_triplets_Sig_Right.prediction],'LineWidth',2);
legend('Scamble','Real');
title('Significant Right-Preferring Triplets');
xlabel('Triplet');
ylabel('Fraction Trials Going Left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,3)
% THIS IS A MISNOMER, GREAT5 IS ACTUALLY GREAT4
diffNum = double([outputTriplets.saveAll_triplets.number_triplets]) - [outputTriplets.saveAll_triplets.shuffle_length_mean];
histogram(diffNum, [round(min(diffNum)):1:round(max(diffNum))],'EdgeColor','none');
title('Real - Shuffle # of Doublets');
xlabel('Triplet');
ylabel('# Trials with Triplet');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,4)
ht1 = histogram([outputTriplets.saveAll_triplets_Sig.prediction], [-0.025:.05:1.025],'EdgeColor','none');
ht1.BinWidth = 0.05;
hold on;
ht2 = histogram([outputTriplets.saveAll_triplets_NotSig.prediction], [-0.025:.05:1.025],'EdgeColor','none');
ht2.BinWidth = 0.05;
xlabel('Fraction Trials Going Left');
title('Significant Triplets - Histogram');
ylabel('# of Triplets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,5)
histogram(sort([outputTriplets.saveAll_triplets_Sig_Left.prediction]-[outputTriplets.saveAll_triplets_Sig_Left.triplet_shuffle]),[-1.025:.05:1.025],'EdgeColor','none');
hold on
histogram(sort([outputTriplets.saveAll_triplets_Sig_Right.prediction]-[outputTriplets.saveAll_triplets_Sig_Right.triplet_shuffle]),[-1.025:.05:1.025],'EdgeColor','none');
xlim([-1 1]);
legend('Left Triplets', 'Right Triplets');
title('Real - Shuffle Prediction');
xlabel('Doublet');
ylabel('Difference Fraction Trials Going Left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,6)
h1 = histogram([outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict], [-0.05:.1:1.05],'EdgeColor','none');
hold on
h2 = histogram([outputTriplets.saveAll_triplets_Sig_Right.prediction], [-0.05:.1:1.05],'EdgeColor','none');
h3 = histogram([outputTriplets.saveAll_triplets_Sig_Right.cell3only_predict],[-0.05:.1:1.05],'EdgeColor','none');
xbins = h1.BinEdges+(h1.BinWidth/2);
xbins = xbins(1:end-1);
plot(xbins,h1.Values,'bo-');
plot(xbins,h2.Values,'ro-');
plot(xbins,h3.Values,'yo-');
legend('1 and 2 only','1, 2, and 3','3 only');
xlabel('Fraction Trials Going Left');
ylabel('# of Sub-Triplets');
title('Triplets - Right Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,9)
h1 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [-0.05:.1:1.05],'EdgeColor','none');
hold on
h2 = histogram([outputTriplets.saveAll_triplets_Sig_Left.prediction], [-0.05:.1:1.05],'EdgeColor','none');
h3 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict],[-0.05:.1:1.05],'EdgeColor','none');
xbins = h1.BinEdges+(h1.BinWidth/2);
xbins = xbins(1:end-1);
plot(xbins,h1.Values,'bo-');
plot(xbins,h2.Values,'ro-');
plot(xbins,h3.Values,'yo-');
legend('1 and 2 only','1, 2, and 3','3 only');
xlabel('Fraction Trials Going Left');
ylabel('# of Sub-Triplets');
title('Triplets - Left Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
% 
% subplot(4,3,7)
% histogram([outputTriplets.saveAll_triplets_Sig_Right(:).cell1and2only_predict],10,'EdgeColor','none');
% hold on;
% histogram([outputTriplets.saveAll_triplets_Sig_Right(:).cell2and3only_predict],10,'EdgeColor','none');
% % ylim([0 1500]);
% legend('Cell 1 and 2 of Triplet ONLY', 'Cell 2 and 3 of Triplet ONLY','Location','northwest');
% title('Triplets - Right Predicting');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% set(gca,'box','off')
% set(gca, 'FontName', 'Arial')
% 
% subplot(4,3,8)
% histogram([outputTriplets.saveAll_triplets_Sig_Left(:).cell1and2only_predict],10,'EdgeColor','none');
% hold on;
% histogram([outputTriplets.saveAll_triplets_Sig_Left(:).cell2and3only_predict],10,'EdgeColor','none');
% % ylim([0 1500]);
% legend('Cell 1 and 2 of Triplet ONLY', 'Cell 2 and 3 of Triplet ONLY','Location','northwest');
% title('Triplets - Left Predicting');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% set(gca,'box','off')
% set(gca, 'FontName', 'Arial')


% subplot(4,3,10)
% h1 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [-0.05:.1:1.05],'EdgeColor','none');
% hold on
% h2 = histogram([outputTriplets.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle], [-0.05:.1:1.05],'EdgeColor','none');
% % h3 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict],[-0.05:.1:1.05],'EdgeColor','none')
% xbins = h1.BinEdges+(h1.BinWidth/2);
% xbins = xbins(1:end-1);
% plot(xbins,h1.Values,'bo-');
% plot(xbins,h2.Values,'ro-');
% % plot(xbins,h3.Values,'yo-');
% legend('1 and 2 only','1 and 2 only - shuffle');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% title('Triplets - Left Predicting');
% set(gca,'box','off')
% set(gca, 'FontName', 'Arial')

% 
% subplot(4,3,11)
% h1 = histogram([outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict], [-0.05:.1:1.05],'EdgeColor','none');
% hold on
% h2 = histogram([outputTriplets.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle], [-0.05:.1:1.05],'EdgeColor','none');
% % h3 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict],[-0.05:.1:1.05],'EdgeColor','none')
% xbins = h1.BinEdges+(h1.BinWidth/2);
% xbins = xbins(1:end-1);
% plot(xbins,h1.Values,'bo-');
% plot(xbins,h2.Values,'ro-');
% % plot(xbins,h3.Values,'yo-');
% legend('1 and 2 only','1 and 2 only - shuffle');
% xlabel('Fraction Trials Going Left');
% ylabel('# of Sub-Triplets');
% title('Triplets - Right Predicting');
% set(gca,'box','off')
% set(gca, 'FontName', 'Arial')


subplot(4,3,10)
scatter([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [outputTriplets.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle], 10,'bo','filled','MarkerFaceAlpha',.3);
hold on;
plot([0 1], [0 1],'r')
xlabel('cell 1 and 2 only - Real');
ylabel('cell 1 and 2 only - Shuffle');
title('Left - % predict left - Real v Shuffle');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,11)
scatter(1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict], 1-[outputTriplets.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle], 10,'ro','filled','MarkerFaceAlpha',.3);
hold on;
plot([0 1], [0 1],'k')
xlabel('cell 1 and 2 only - Real');
ylabel('cell 1 and 2 only - Shuffle');
title('Right - % predict right - Real v Shuffle');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')



subplot(4,3,7)
Left_1_2_3     = [outputTriplets.saveAll_triplets_Sig_Left.prediction];
Left_1_2_no3   = [outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict];
Left_no1_no2_3 = [outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict];
outputBarSEM_Left = nieh_barSEM(Left_1_2_3, Left_1_2_no3, Left_no1_no2_3)
xticklabels({'1,2,3','1,2,no3','3 only'});
ylabel('Fraction Trials Going Left');
title('Triplets - Left Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,8)
Right_1_2_3     = 1-[outputTriplets.saveAll_triplets_Sig_Right.prediction];
Right_1_2_no3   = 1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict];
Right_no1_no2_3 = 1-[outputTriplets.saveAll_triplets_Sig_Right.cell3only_predict];
outputBarSEM_Right = nieh_barSEM(Right_1_2_3, Right_1_2_no3, Right_no1_no2_3)
xticklabels({'1,2,3','1,2,no3','3 only'});
ylabel('Fraction Trials Going Right');


title('Triplets - Right Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


%%
figure;
subplot(1,2,1)
Left_1_2_3     = [outputTriplets.saveAll_triplets_Sig_Left.prediction];
Left_1_2_no3   = [outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict];
Left_no1_no2_3 = [outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict];

h1=boxplot([Left_1_2_3; Left_1_2_no3; Left_no1_no2_3]','Whisker',Inf,'PlotStyle','traditional','Widths',.3,'Colors',[.5 .7 1]);
hold on;

set(h1,{'linew'},{2})

scatter(ones(size(Left_1_2_3)).*(1+(rand(size(Left_1_2_3))-0.5)/5),Left_1_2_3,10,'k','filled','MarkerFaceAlpha',0.2);
scatter(ones(size(Left_1_2_no3)).*(2+(rand(size(Left_1_2_no3))-0.5)/5),Left_1_2_no3,10,'k','filled','MarkerFaceAlpha',0.2);
scatter(ones(size(Left_no1_no2_3)).*(3+(rand(size(Left_no1_no2_3))-0.5)/5),Left_no1_no2_3,10,'k','filled','MarkerFaceAlpha',0.2);

xticklabels({'1,2,3','1,2,no3','3 only'});
ylabel('Fraction Trials Going Left');
title('Triplets - Left Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')





subplot(1,2,2)
Right_1_2_3     = 1-[outputTriplets.saveAll_triplets_Sig_Right.prediction];
Right_1_2_no3   = 1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict];
Right_no1_no2_3 = 1-[outputTriplets.saveAll_triplets_Sig_Right.cell3only_predict];


h1=boxplot([Right_1_2_3; Right_1_2_no3; Right_no1_no2_3]','Whisker',Inf,'PlotStyle','traditional','Widths',.3,'Colors','r');
hold on;

set(h1,{'linew'},{2})

scatter(ones(size(Right_1_2_3)).*(1+(rand(size(Right_1_2_3))-0.5)/5),Right_1_2_3,10,'k','filled','MarkerFaceAlpha',0.2);
scatter(ones(size(Right_1_2_no3)).*(2+(rand(size(Right_1_2_no3))-0.5)/5),Right_1_2_no3,10,'k','filled','MarkerFaceAlpha',0.2);
scatter(ones(size(Right_no1_no2_3)).*(3+(rand(size(Right_no1_no2_3))-0.5)/5),Right_no1_no2_3,10,'k','filled','MarkerFaceAlpha',0.2);



xticklabels({'1,2,3','1,2,no3','3 only'});
ylabel('Fraction Trials Going Right');
title('Triplets - Right Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


%% Violin Plots

figure;
subplot(1,2,1)
Left_1_2_3     = [outputTriplets.saveAll_triplets_Sig_Left.prediction];
Left_1_2_no3   = [outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict];
Left_no1_no2_3 = [outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict];

h1=distributionPlot([Left_1_2_3; Left_1_2_no3; Left_no1_no2_3]','showMM',0,'globalNorm',1,'histOpt',0,'divFactor',[0:.1:1]);
hold on;
h1b=boxplot([Left_1_2_3; Left_1_2_no3; Left_no1_no2_3]','Whisker',Inf,'PlotStyle','traditional','Widths',.08,'Colors',[.5 .7 1]);
set(h1b,{'linew'},{2})

xticklabels({'1,2,3','1,2,no3','3 only'});
ylabel('Fraction Trials Going Left');
title('Triplets - Left Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
%axis square

subplot(1,2,2)
Right_1_2_3     = 1-[outputTriplets.saveAll_triplets_Sig_Right.prediction];
Right_1_2_no3   = 1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict];
Right_no1_no2_3 = 1-[outputTriplets.saveAll_triplets_Sig_Right.cell3only_predict];


h2=distributionPlot([Right_1_2_3; Right_1_2_no3; Right_no1_no2_3]','showMM',0,'globalNorm',1,'histOpt',0,'divFactor',[0:.1:1]);
hold on;
h2b=boxplot([Right_1_2_3; Right_1_2_no3; Right_no1_no2_3]','Whisker',Inf,'PlotStyle','traditional','Widths',.08,'Colors','r');
set(h2b,{'linew'},{2})

xticklabels({'1,2,3','1,2,no3','3 only'});
ylabel('Fraction Trials Going Right');
title('Triplets - Right Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
%axis square
