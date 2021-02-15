function plotDoublets(outputDoublets)

rng(1)

%% Plot the data in Fig. 4

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
direction_real_shuffle = [outputDoublets.saveAll_doublets.asymmetry_length_diff; outputDoublets.saveAll_doublets.mean_asymmetry_diff_length_shuf];
diff_sym = [direction_real_shuffle(1,:)-direction_real_shuffle(2,:)];
[~,p_diff_sym] = ttest(direction_real_shuffle(1,:), direction_real_shuffle(2,:));
disp(['paired t test : p=' num2str(p_diff_sym)]); 
h=boxplot(diff_sym,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
set(h,{'linew'},{2})
hold on
scatter(ones(size(diff_sym)).*(1+(rand(size(diff_sym))-0.5)/5),diff_sym,5,'k','filled','MarkerFaceAlpha',0.2);
ylim([-10 28])
xticklabels({'Real-Shuffle'});
ylabel('# Trials Doublet Appears (Forwards-Backwards)');
title(['Fig. 4c, n=' num2str(length(diff_sym)) ' doublets, p=' num2str(signrank(diff_sym))]);
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,3)
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


%% Plot the data in Extended Data Fig. 5

figure;

subplot(2,2,1)
diffNum = [outputDoublets.saveAll_doublets.number_doublets] - [outputDoublets.saveAll_doublets.shuffle_length];
histogram(diffNum, [round(min(diffNum)):1:30],'EdgeColor','none');
title('Real - Shuffle # of Doublets');
xlabel('# trials doublet appears');
ylabel('# of doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,2)
diffAsym = [outputDoublets.saveAll_doublets.asymmetry_length_diff] - [outputDoublets.saveAll_doublets.mean_asymmetry_diff_length_shuf];
histogram(diffAsym, [round(min(diffAsym)):1:30],'EdgeColor','none');
title('Asymmetry of Doublets (Real - Shuf of # Forward - # Backwards)');
xlabel('Directionality index (real - shuf)');
ylabel('# of doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,3)
histogram([outputDoublets.saveAll_doublets_Sig_Left.prediction]-[outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle],[-1.025:.05:1.025],'EdgeColor','none');
title('Left Doublets, Real - Shuffle Prediction');
xlabel('Fraction doublets went left');
ylabel('# of Doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,4)
histogram(-1*([outputDoublets.saveAll_doublets_Sig_Right.prediction]-[outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle]),[-1.025:.05:1.025],'EdgeColor','none','FaceColor','r');
title('Right Doublets, Real - Shuffle Prediction');
xlabel('Fraction doublets went right');
ylabel('# of Doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


