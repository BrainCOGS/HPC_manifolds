function findDoublets_20190509_plotter(outputDoublets)

rng(1)

%% Plot the data
figure('units','normalized','outerposition',[0 0 1 1])


subplot(4,3,1)
plot([outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle],'LineWidth',2);
hold on;
plot([outputDoublets.saveAll_doublets_Sig_Left.prediction],'LineWidth',2);
plot([outputDoublets.saveAll_doublets_Sig_Left.only_first_predict],'LineWidth',2);
plot([outputDoublets.saveAll_doublets_Sig_Left.only_second_predict],'LineWidth',2);
legend('Scamble','Real','Cell1','Cell2');
title('Significant Left-Preferring Doublets');
xlabel('Doublet');
ylabel('Fraction Trials Going Left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,2)
plot([outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle],'LineWidth',2);
hold on;
plot([outputDoublets.saveAll_doublets_Sig_Right.prediction],'LineWidth',2);
plot([outputDoublets.saveAll_doublets_Sig_Right.only_first_predict],'LineWidth',2);
plot([outputDoublets.saveAll_doublets_Sig_Right.only_second_predict],'LineWidth',2);
legend('Scamble','Real','Cell1','Cell2');
%legend('Scamble','Real');
title('Significant Right-Preferring Doublets');
xlabel('Doublet');
ylabel('Fraction Trials Going Left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,3)
histogram(sort([outputDoublets.saveAll_doublets_Sig_Left.prediction]-[outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle]),[-1.025:.05:1.025],'EdgeColor','none');
hold on
histogram(sort([outputDoublets.saveAll_doublets_Sig_Right.prediction]-[outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle]),[-1.025:.05:1.025],'EdgeColor','none');
legend('Left Doublets', 'Right Doublets');
title('Real - Shuffle Prediction');
xlabel('Doublet % went Left');
ylabel('Difference Fraction Trials Going Left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

% subplot(2,3,3)
% plot(sort([saveAll_doublets_Sig_Left.prediction]-[saveAll_doublets_Sig_Left.doublet_shuffle]),'LineWidth',2);
% hold on
% plot(sort([saveAll_doublets_Sig_Right.prediction]-[saveAll_doublets_Sig_Right.doublet_shuffle]),'LineWidth',2);
% legend('Left Doublets', 'Right Doublets');
% title('Real - Shuffle Prediction');
% xlabel('Doublet');
% ylabel('Difference Fraction Trials Going Left');

subplot(4,3,4)
diffNum = [outputDoublets.saveAll_doublets.number_doublets] - [outputDoublets.saveAll_doublets.shuffle_length];
histogram(diffNum, [round(min(diffNum)):1:30],'EdgeColor','none');
%histogram(diffNum, [round(min(diffNum)):2:round(max(diffNum))],'EdgeColor','none');
hold on;
numDoublets = [outputDoublets.saveAll_doublets_Sig.number_doublets];
histogram(numDoublets, [round(min(numDoublets)):1:30],'EdgeColor','none');
%histogram(numDoublets, [round(min(numDoublets)):1:round(max(numDoublets))],'EdgeColor','none');
title('Real - Shuffle # of Doublets');
xlabel('Doublet');
ylabel('# Trials with Doublet');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,5)
h1 = histogram([outputDoublets.saveAll_doublets_Sig.prediction],[-0.025:.05:1.025],'EdgeColor','none');
hold on
h2 = histogram([outputDoublets.saveAll_doublets_NotSig.prediction],[-0.025:.05:1.025],'EdgeColor','none');
xlabel('Fraction Trials Going Left');
title('Significant Doublets - Histogram');
ylabel('# of Doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,6)
%diffAsym = [outputDoublets.saveAll_doublets.number_doublets] - [outputDoublets.saveAll_doublets.asymmetry_length];
diffAsym = [outputDoublets.saveAll_doublets.asymmetry_length_diff] - [outputDoublets.saveAll_doublets.mean_asymmetry_diff_length_shuf];
histogram(diffAsym, [round(min(diffAsym)):1:30],'EdgeColor','none');
%histogram(diffAsym, [round(min(diffAsym)):2:round(max(diffAsym))],'EdgeColor','none');
%title('Asymmetry of Doublets (# Forward - # Backwards)');
title('Asymmetry of Doublets (Real - Shuf of # Forward - # Backwards)');
xlabel('Doublet');
ylabel('Number of Trials with Doublet');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(4,3,9)
scatter([outputDoublets.saveAll_doublets.number_doublets], [outputDoublets.saveAll_doublets.asymmetry_length],10,'bo','filled','MarkerFaceAlpha',.3);
title('Asymmetry of Doublets');
xlabel('Forwards')
ylabel('Backwards')
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
axis square
xlim([0 40])
ylim([0 40])
hold on
plot([0 40],[0 40],'r');

subplot(4,3,7)
scatter([outputDoublets.saveAll_doublets_Sig_Left.prediction], [outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle],10,'bo','filled','MarkerFaceAlpha',.3);
hold on
title('Real vs Shuffle Prediction - Left');
xlabel('Doublet % went Left');
ylabel('Shuffle Doublet % went left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
plot([0 1],[0 1],'k');
axis square

subplot(4,3,8)
scatter(1-[outputDoublets.saveAll_doublets_Sig_Right.prediction], 1-[outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle],10,'ro','filled','MarkerFaceAlpha',.3);
hold on
title('Real vs Shuffle Prediction - Right');
xlabel('Doublet % went Left');
ylabel('Shuffle Doublet % went left');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
plot([0 1],[0 1],'k');
axis square


subplot(4,3,10)
scatter([outputDoublets.saveAll_doublets.number_doublets],[outputDoublets.saveAll_doublets.shuffle_length],10,'ro','filled','MarkerFaceAlpha',.3);
hold on
title('Real vs Shuffle Prediction - Length');
xlabel('# Real Doublets');
ylabel('# Shuffle Doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
plot([0 40],[0 40],'k');
axis square
% 
% backwards = [outputDoublets.saveAll_doublets_Sig.asymmetry_length];
% forward = [outputDoublets.saveAll_doublets_Sig.number_doublets];
% figure; histogram2(forward, backwards,[0:1:20],[0:1:20],'DisplayStyle','tile','ShowEmptyBins','on', 'EdgeColor','none'); 
% axis square;
% % colormap(flipud(gray))
% 
% realnum = [outputDoublets.saveAll_doublets_great5.number_doublets];
% shufnum = [outputDoublets.saveAll_doublets_great5.shuffle_length];
% figure; histogram2(realnum, shufnum,[0:1:20],[0:1:20],'DisplayStyle','tile','ShowEmptyBins','on', 'EdgeColor','none'); 
% axis square;

subplot(4,3,11)
diffNum = [outputDoublets.saveAll_doublets.number_doublets] - [outputDoublets.saveAll_doublets.shuffle_length];
diffAsym = [outputDoublets.saveAll_doublets.number_doublets] - [outputDoublets.saveAll_doublets.asymmetry_length];
nieh_barSEM(diffNum, diffAsym);
legend('Difference # Doublets','Difference - Assymmetry');
ylim([0 7]);

subplot(4,3,12)
diffLeft = [outputDoublets.saveAll_doublets_Sig_Left.prediction]-[outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle];
diffRight = (1-[outputDoublets.saveAll_doublets_Sig_Right.prediction])-(1-[outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle]);
nieh_barSEM(diffLeft, diffRight);



%% Some new stuff from second draft that PIs saw


% 
% figure; 
% 
% numberDub_real_shuffle = [outputDoublets.saveAll_doublets.number_doublets; outputDoublets.saveAll_doublets.shuffle_length];
% subplot(2,2,1)
% p1 = plot(numberDub_real_shuffle,'o-', 'Color',[0 0 0 .1]);
% xlim([0 3])
% xticks([1 2])
% xticklabels({'Real' ,'Shuffle'});
% ylabel('# Trials Doublet Appears');
% title('Same as Fig. 4b, n=16088 doublets');
% 
% direction_real_shuffle = [outputDoublets.saveAll_doublets.number_doublets; outputDoublets.saveAll_doublets.asymmetry_length];
% subplot(2,2,2)
% plot(direction_real_shuffle,'o-', 'Color',[0 0 0 .1])
% xlim([0 3])
% xticks([1 2])
% xticklabels({'Forwards' ,'Backwards'});
% ylabel('# Trials Doublet Appears');
% title('Same as Fig. 4c, n=16088 doublets');
% 
% 
% left_real_shuffle = [outputDoublets.saveAll_doublets_Sig_Left.prediction; outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle];
% subplot(2,2,3)
% plot(left_real_shuffle,'o-','Color',[0 0 1 .1])
% xlim([0 3])
% xticks([1 2])
% xticklabels({'Real' ,'Shuffle'});
% ylabel('% Trials Went Left');
% title('Same as Fig. 4f (left) n=922 doublets');
% 
% right_real_shuffle = 1-[outputDoublets.saveAll_doublets_Sig_Right.prediction; outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle];
% subplot(2,2,4)
% plot(right_real_shuffle,'o-', 'Color',[1 0 0 .1])
% xlim([0 3])
% xticks([1 2])
% xticklabels({'Real' ,'Shuffle'});
% ylabel('% Trials Went Right');
% title('Same as Fig. 4f (right) n=1227 doublets');
% 


figure; 
subplot(2,2,1)
numberDub_real_shuffle = [outputDoublets.saveAll_doublets.number_doublets; outputDoublets.saveAll_doublets.shuffle_length];
diff_num = [numberDub_real_shuffle(1,:)-numberDub_real_shuffle(2,:)];
[~,p_diff_num] = ttest(numberDub_real_shuffle(1,:), numberDub_real_shuffle(2,:));
disp(['paired t test : p=' num2str(p_diff_num)]); 
h=boxplot(diff_num,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
set(h,{'linew'},{2})
hold on
scatter(ones(size(diff_num)).*(1+(rand(size(diff_num))-0.5)/5),diff_num,5,'k','filled','MarkerFaceAlpha',0.2);
ylim([-10 25])
xticklabels({'Real-Shuffle'});
ylabel('# Trials Doublet Appears');
title(['Fig. 4b, n=' num2str(length(diff_num)) ' doublets, p=' num2str(signrank(diff_num))]);
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,2)
%direction_real_shuffle = [outputDoublets.saveAll_doublets.number_doublets; outputDoublets.saveAll_doublets.asymmetry_length];
direction_real_shuffle = [outputDoublets.saveAll_doublets.asymmetry_length_diff; outputDoublets.saveAll_doublets.mean_asymmetry_diff_length_shuf];
diff_sym = [direction_real_shuffle(1,:)-direction_real_shuffle(2,:)];
[~,p_diff_sym] = ttest(direction_real_shuffle(1,:), direction_real_shuffle(2,:));
disp(['paired t test : p=' num2str(p_diff_sym)]); 
h=boxplot(diff_sym,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
set(h,{'linew'},{2})
hold on
scatter(ones(size(diff_sym)).*(1+(rand(size(diff_sym))-0.5)/5),diff_sym,5,'k','filled','MarkerFaceAlpha',0.2);
ylim([-10 28])
%xticklabels({'Forwards-Backwards'});
xticklabels({'Real-Shuffle'});
ylabel('# Trials Doublet Appears (Forwards-Backwards)');
title(['Fig. 4c, n=' num2str(length(diff_sym)) ' doublets, p=' num2str(signrank(diff_sym))]);
set(gca,'box','off')
set(gca, 'FontName', 'Arial')



subplot(2,2,3)
%plotSpread({left_real_shuffle(1,:)-left_real_shuffle(2,:)},'distributionMarkers',10,'SpreadWidth',.2)
left_real_shuffle = [outputDoublets.saveAll_doublets_Sig_Left.prediction; outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle];
diff_left = [left_real_shuffle(1,:)-left_real_shuffle(2,:)];
[~,p_diff_left] = ttest(left_real_shuffle(1,:), left_real_shuffle(2,:));
disp(['paired t test : p=' num2str(p_diff_left)]); 
h=boxplot(diff_left,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
set(h,{'linew'},{2})
hold on
scatter(ones(size(diff_left)).*(1+(rand(size(diff_left))-0.5)/5),diff_left,5,'k','filled','MarkerFaceAlpha',0.2);
ylim([-0.2 .7])
xticklabels({'Real-Shuffle'});
ylabel('% Trials Went Left');
title(['Fig. 4f (left) n=' num2str(length(diff_left)) ' doublets, p=' num2str(p_diff_left)]);
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,4)
right_real_shuffle = 1-[outputDoublets.saveAll_doublets_Sig_Right.prediction; outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle];
diff_right = [right_real_shuffle(1,:)-right_real_shuffle(2,:)];
[~,p_diff_right] = ttest(right_real_shuffle(1,:), right_real_shuffle(2,:));
disp(['paired t test : p=' num2str(p_diff_right)]); 
h=boxplot(diff_right,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
set(h,{'linew'},{2})
hold on
scatter(ones(size(diff_right)).*(1+(rand(size(diff_right))-0.5)/5),diff_right,5,'k','filled','MarkerFaceAlpha',0.2);
ylim([-0.2 .7])
xticklabels({'Real-Shuffle'});
ylabel('% Trials Went Right');
title(['Fig. 4f (right) n=' num2str(length(diff_right)) ' doublets, p=' num2str(p_diff_right)]);
set(gca,'box','off')
set(gca, 'FontName', 'Arial')












