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

leftSig_triplet(:,2) = Left_1_2_no3;
leftSig_triplet(:,1) = Left_1_2_3; 
leftSig_triplet(:,3) = Left_no1_no2_3;
[~,pl1] = ttest(leftSig_triplet(:,1), leftSig_triplet(:,2));
[~,pl2] = ttest(leftSig_triplet(:,1), leftSig_triplet(:,3));
[~,pl3] = ttest(leftSig_triplet(:,2), leftSig_triplet(:,3));
disp(['Fig. S9g col 1 vs col 2, paired t test : p=' num2str(pl1*3) ' n=' num2str(size(leftSig_triplet,1))]);
disp(['Fig. S9g col 1 vs col 3, paired t test : p=' num2str(pl2*3) ' n=' num2str(size(leftSig_triplet,1))]);
disp(['Fig. S9g col 2 vs col 3, paired t test : p=' num2str(pl3*3) ' n=' num2str(size(leftSig_triplet,1))]);

subplot(2,2,2)
scatter([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [outputTriplets.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle], 10,'bo','filled','MarkerFaceAlpha',.3);
hold on;
plot([0 1], [0 1],'r')
xlabel('cell 1 and 2 only - Real');
ylabel('cell 1 and 2 only - Shuffle');
axis square
title('Left - fraction predict left - Real v Shuffle');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

[~, pl] = ttest([outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict], [outputTriplets.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle]);
disp(['Fig. S9h Real v Shuf, paired t test : p=' num2str(pl) ' n=' num2str(size(leftSig_triplet,1))]);

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

rightSig_triplet(:,2) = Right_1_2_no3;
rightSig_triplet(:,1) = Right_1_2_3; 
rightSig_triplet(:,3) = Right_no1_no2_3;
[~,pr1] = ttest(rightSig_triplet(:,1), rightSig_triplet(:,2));
[~,pr2] = ttest(rightSig_triplet(:,1), rightSig_triplet(:,3));
[~,pr3] = ttest(rightSig_triplet(:,2), rightSig_triplet(:,3));
disp(['Fig. S9i col 1 vs col 2, paired t test : p=' num2str(pr1*3) ' n=' num2str(size(rightSig_triplet,1))]);
disp(['Fig. S9i col 1 vs col 3, paired t test : p=' num2str(pr2*3) ' n=' num2str(size(rightSig_triplet,1))]);
disp(['Fig. S9i col 2 vs col 3, paired t test : p=' num2str(pr3*3) ' n=' num2str(size(rightSig_triplet,1))]);

subplot(2,2,4)
scatter(1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict], 1-[outputTriplets.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle], 10,'ro','filled','MarkerFaceAlpha',.3);
hold on;
plot([0 1], [0 1],'k')
xlabel('cell 1 and 2 only - Real');
ylabel('cell 1 and 2 only - Shuffle');
axis square
title('Right - fraction predict right - Real v Shuffle');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

[~, pr] = ttest([outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict], [outputTriplets.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle]);
disp(['Fig. S9j Real v Shuf, paired t test : p=' num2str(pr) ' n=' num2str(size(rightSig_triplet,1))]);

figure;
subplot(1,2,1)
diffLeft = [outputTriplets.saveAll_triplets_Sig_Left.cell1and2only_predict] - [outputTriplets.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle];
scatter(ones(size(diffLeft)).*(1+(rand(size(diffLeft))-0.5)/5),diffLeft,5,'k','filled','MarkerFaceAlpha',0.2);
hold on;
h=boxplot(diffLeft,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','b');
set(h,{'linew'},{2})
xticklabels({'Real - Shuffle'});
ylabel('Fraction predict left');
title('Left triplets')
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(1,2,2)
diffRight = (1-[outputTriplets.saveAll_triplets_Sig_Right.cell1and2only_predict]) - (1-[outputTriplets.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle]);
scatter(ones(size(diffRight)).*(1+(rand(size(diffRight))-0.5)/5),diffRight,5,'k','filled','MarkerFaceAlpha',0.2);
hold on;
h=boxplot(diffRight,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','b');
set(h,{'linew'},{2})
xticklabels({'Real - Shuffle'});
ylabel('Fraction predict right');
title('Right triplets')
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


%%

figure;

subplot(3,2,1)
diffNum = double([outputTriplets.saveAll_triplets.number_triplets]) - [outputTriplets.saveAll_triplets.shuffle_length_mean];
histogram(diffNum, [round(min(diffNum)):1:round(max(diffNum))],'EdgeColor','none');
title('Real - Shuffle # of Triplets');
xlabel('# trials where triplet appears (real - shuffle)');
ylabel('# of triplets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,2)
scatter(ones(size(diffNum)).*(1+(rand(size(diffNum))-0.5)/5),diffNum,5,'k','filled','MarkerFaceAlpha',0.2);
hold on;
h=boxplot(diffNum,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','b');
set(h,{'linew'},{2})
xticklabels({'Real - Shuffle'});
ylabel('# trials triplet appears');
set(gca,'box','off')

[~, p_diffNum] = ttest(diffNum);
disp(['Fig. S9e Real v Shuf, paired t test : p=' num2str(p_diffNum) ' n=' num2str(size(diffNum,2))]);

subplot(3,2,3)
diffPred_left = [outputTriplets.saveAll_triplets_Sig_Left.prediction]-[outputTriplets.saveAll_triplets_Sig_Left.triplet_shuffle];
histogram(diffPred_left,[-1.025:.05:1.025],'EdgeColor','none');
xlim([-.5 1]);
title('Left Triplets, Real - Shuffle Prediction');
xlabel('Fraction trials went left');
ylabel('# of Triplets');
ylim([0 600])
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,4)
scatter(ones(size(diffPred_left)).*(1+(rand(size(diffPred_left))-0.5)/5),diffPred_left,5,'k','filled','MarkerFaceAlpha',0.2);
hold on;
h=boxplot(diffPred_left,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','b');
set(h,{'linew'},{2})
xticklabels({'Real - Shuffle'});
ylabel('Fraction trials went left');
set(gca,'box','off')
ylim([-.5 1]);

[~, p_diffPred_left] = ttest(diffPred_left);
disp(['Fig. S9f (left) Real v Shuf, paired t test : p=' num2str(p_diffPred_left) ' n=' num2str(size(diffPred_left,2))]);

subplot(3,2,5)
diffPred_right = -1*([outputTriplets.saveAll_triplets_Sig_Right.prediction]-[outputTriplets.saveAll_triplets_Sig_Right.triplet_shuffle]);
histogram(diffPred_right,[-1.025:.05:1.025],'EdgeColor','none');
xlim([-.5 1]);
title('Right Triplets, Real - Shuffle Prediction');
xlabel('Fraction trials went right');
ylabel('# of Triplets');
ylim([0 600])
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(3,2,6)
scatter(ones(size(diffPred_right)).*(1+(rand(size(diffPred_right))-0.5)/5),diffPred_right,5,'k','filled','MarkerFaceAlpha',0.2);
hold on;
h=boxplot(diffPred_right,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','b');
set(h,{'linew'},{2})
xticklabels({'Real - Shuffle'});
ylabel('Fraction trials went right');
set(gca,'box','off')
ylim([-.5 1]);

[~, p_diffPred_right] = ttest(diffPred_right);
disp(['Fig. S9f (right) Real v Shuf, paired t test : p=' num2str(p_diffPred_right) ' n=' num2str(size(diffPred_right,2))]);

set(gcf, 'Position', [500, 90, 280, 680])


