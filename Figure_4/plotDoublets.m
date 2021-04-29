function outputPlotDoublets = plotDoublets(outputDoublets)

rng(1)

%% Make upfront calculations

numberDub_real_shuffle = [outputDoublets.saveAll_doublets.number_doublets; outputDoublets.saveAll_doublets.shuffle_length];
diff_num = [numberDub_real_shuffle(1,:)-numberDub_real_shuffle(2,:)];
[~,p_diff_num] = ttest(numberDub_real_shuffle(1,:), numberDub_real_shuffle(2,:));
disp(['Fig. 4b, paired t test : p=' num2str(p_diff_num)]); 

direction_real_shuffle = [outputDoublets.saveAll_doublets.asymmetry_length_diff; outputDoublets.saveAll_doublets.mean_asymmetry_diff_length_shuf];
diff_sym = [direction_real_shuffle(1,:)-direction_real_shuffle(2,:)];
[~,p_diff_sym] = ttest(direction_real_shuffle(1,:), direction_real_shuffle(2,:));
disp(['Fig. 4c, paired t test : p=' num2str(p_diff_sym)]); 

left_real_shuffle = [outputDoublets.saveAll_doublets_Sig_Left.prediction; outputDoublets.saveAll_doublets_Sig_Left.doublet_shuffle];
diff_left = [left_real_shuffle(1,:)-left_real_shuffle(2,:)];
[~,p_diff_left] = ttest(left_real_shuffle(1,:), left_real_shuffle(2,:));
disp(['Fig 4f (left), paired t test : p=' num2str(p_diff_left)]); 

right_real_shuffle = 1-[outputDoublets.saveAll_doublets_Sig_Right.prediction; outputDoublets.saveAll_doublets_Sig_Right.doublet_shuffle];
diff_right = [right_real_shuffle(1,:)-right_real_shuffle(2,:)];
[~,p_diff_right] = ttest(right_real_shuffle(1,:), right_real_shuffle(2,:));
disp(['Fig 4f (right), paired t test : p=' num2str(p_diff_right)]); 



%% Plot the data in Fig. 4

figure; 

subplot(2,2,1)
h=boxplot(diff_num,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
sourceData_4b = diff_num';
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
h=boxplot(diff_sym,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
sourceData_4c = diff_sym';
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
h=boxplot(diff_left,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
sourceData_4f_left = diff_left'*100;
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
h=boxplot(diff_right,'Whisker',1000,'PlotStyle','traditional','Widths',.3,'Colors','k');
sourceData_4f_right = diff_right'*100;
set(h,{'linew'},{2})
hold on
scatter(ones(size(diff_right)).*(1+(rand(size(diff_right))-0.5)/5),diff_right,5,'k','filled','MarkerFaceAlpha',0.2);
ylim([-0.2 .7])
xticklabels({'Real-Shuffle'});
ylabel('% Trials Went Right');
title(['Fig. 4f (right) n=' num2str(length(diff_right)) ' doublets, p=' num2str(p_diff_right)]);
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


%% Plot the data in Extended Data Fig. 9

figure;

subplot(2,2,1)
h1 = histogram(diff_num, [round(min(diff_num)):1:30],'EdgeColor','none');
sourceData_S9a_x = h1.BinEdges';
sourceData_S9a_y = h1.BinCounts';
title('Real - Shuffle # of Doublets');
xlabel('# trials doublet appears');
ylabel('# of doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,2)
h2 = histogram(diff_sym, [round(min(diff_sym)):1:30],'EdgeColor','none');
sourceData_S9b_x = h2.BinEdges';
sourceData_S9b_y = h2.BinCounts';
title('Asymmetry of Doublets (Real - Shuf of # Forward - # Backwards)');
xlabel('Directionality index (real - shuf)');
ylabel('# of doublets');
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,3)
h3 = histogram(diff_left,[-0.525:.05:1.025],'EdgeColor','none');
sourceData_S9c_left_x = h3.BinEdges'*100;
sourceData_S9c_left_y = h3.BinCounts';
title('Left Doublets, Real - Shuffle Prediction');
xlabel('Fraction doublets went left');
ylabel('# of Doublets');
xlim([-.5 1])
set(gca,'box','off')
set(gca, 'FontName', 'Arial')

subplot(2,2,4)
h4 = histogram(diff_right,[-0.525:.05:1.025],'EdgeColor','none','FaceColor','r');
sourceData_S9c_right_x = h4.BinEdges'*100;
sourceData_S9c_right_y = h4.BinCounts';
title('Right Doublets, Real - Shuffle Prediction');
xlabel('Fraction doublets went right');
ylabel('# of Doublets');
xlim([-.5 1])
set(gca,'box','off')
set(gca, 'FontName', 'Arial')


%% Save data for output

outputPlotDoublets.sourceData_4b = sourceData_4b;
outputPlotDoublets.sourceData_4c = sourceData_4c;
outputPlotDoublets.sourceData_4f_left = sourceData_4f_left;
outputPlotDoublets.sourceData_4f_right = sourceData_4f_right;
outputPlotDoublets.sourceData_S9a_x = sourceData_S9a_x;
outputPlotDoublets.sourceData_S9a_y = sourceData_S9a_y;
outputPlotDoublets.sourceData_S9b_x = sourceData_S9b_x;
outputPlotDoublets.sourceData_S9b_y = sourceData_S9b_y;
outputPlotDoublets.sourceData_S9c_left_x = sourceData_S9c_left_x;
outputPlotDoublets.sourceData_S9c_left_y = sourceData_S9c_left_y;
outputPlotDoublets.sourceData_S9c_right_x = sourceData_S9c_right_x;
outputPlotDoublets.sourceData_S9c_right_y = sourceData_S9c_right_y;