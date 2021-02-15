function plotTriplets(outputTriplets)

rng(1)

%% Plot data in Fig. 4

figure;

subplot(2,2,1)
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
axis square
title('Triplets - Left Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,2)
scatter([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [outputTriplets.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle], 10,'bo','filled','MarkerFaceAlpha',.3);
hold on;
plot([0 1], [0 1],'r')
xlabel('cell 1 and 2 only - Real');
ylabel('cell 1 and 2 only - Shuffle');
axis square
title('Left - % predict left - Real v Shuffle');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,3)
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
axis square
title('Triplets - Right Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,4)
scatter(1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict], 1-[outputTriplets.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle], 10,'ro','filled','MarkerFaceAlpha',.3);
hold on;
plot([0 1], [0 1],'k')
xlabel('cell 1 and 2 only - Real');
ylabel('cell 1 and 2 only - Shuffle');
axis square
title('Right - % predict right - Real v Shuffle');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


%%

figure;

subplot(3,2,1)
diffNum = double([outputTriplets.saveAll_triplets.number_triplets]) - [outputTriplets.saveAll_triplets.shuffle_length_mean];
histogram(diffNum, [round(min(diffNum)):1:round(max(diffNum))],'EdgeColor','none');
title('Real - Shuffle # of Doublets');
xlabel('# trials where triplet appears (real - shuffle)');
ylabel('# of triplets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,3)
histogram([outputTriplets.saveAll_triplets_Sig_Left.prediction]-[outputTriplets.saveAll_triplets_Sig_Left.triplet_shuffle],[-1.025:.05:1.025],'EdgeColor','none');
xlim([-1 1]);
title('Left Triplets, Real - Shuffle Prediction');
xlabel('Fraction trials went left');
ylabel('# of Triplets');
ylim([0 600])
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,4)
histogram(-1*([outputTriplets.saveAll_triplets_Sig_Right.prediction]-[outputTriplets.saveAll_triplets_Sig_Right.triplet_shuffle]),[-1.025:.05:1.025],'EdgeColor','none');
xlim([-1 1]);
title('Right Triplets, Real - Shuffle Prediction');
xlabel('Fraction trials went right');
ylabel('# of Triplets');
ylim([0 600])
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,5)
h1 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [-0.05:.1:1.05],'EdgeColor','none');
hold on
h2 = histogram([outputTriplets.saveAll_triplets_Sig_Left.prediction], [-0.05:.1:1.05],'EdgeColor','none');
h3 = histogram([outputTriplets.saveAll_triplets_Sig_Left.cell3only_predict],[-0.05:.1:1.05],'EdgeColor','none');
xbins = h1.BinEdges+(h1.BinWidth/2);
xbins = xbins(1:end-1);
plot(xbins,h1.Values,'b-');
plot(xbins,h2.Values,'r-');
plot(xbins,h3.Values,'y-');
legend('1 and 2 only','1, 2, and 3','3 only');
xlabel('Fraction trials went left');
ylabel('# of Triplets');
title('Triplets - Left Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,6)
h1 = histogram(1-([outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict]), [-0.05:.1:1.05],'EdgeColor','none');
hold on
h2 = histogram(1-([outputTriplets.saveAll_triplets_Sig_Right.prediction]), [-0.05:.1:1.05],'EdgeColor','none');
h3 = histogram(1-([outputTriplets.saveAll_triplets_Sig_Right.cell3only_predict]),[-0.05:.1:1.05],'EdgeColor','none');
xbins = h1.BinEdges+(h1.BinWidth/2);
xbins = xbins(1:end-1);
plot(xbins,h1.Values,'b-');
plot(xbins,h2.Values,'r-');
plot(xbins,h3.Values,'y-');
legend('1 and 2 only','1, 2, and 3','3 only');
xlabel('Fraction trials went right');
ylabel('# of Triplets');
title('Triplets - Right Predicting');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


