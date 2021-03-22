function plotSequenceTrialROI_multiple(outputDoublets, fname, out_Y, randDoublet)

nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [0 0], fname,'none','towers',1,1);
trialn = nic_output.behavioralVariables.Trial;
[data_trial_nrs,~,~] = unique(trialn);

%% Make the plots

figure;
hold on;

for k=1:length(randDoublet)
    
    curDoublet = randDoublet(k);
    ax{k} = subplot(5,5,k);
    roi1 = outputDoublets.saveAll_doublets(curDoublet).cells;
    roi1_all(k,:) = roi1;
    numtrial = outputDoublets.saveAll_doublets(curDoublet).trials_appear;
    nic_output = extractVariables(roi1, 2, data_trial_nrs(numtrial)', 2, 100, 0, 5, [0 0], fname,'none','towers',1,1);
    ROI_trial = reshape(nic_output.ROIactivities,[101 length(numtrial) length(roi1)]);

    for j=1:length(numtrial)
        for i=1:length(roi1)
            plotdataY(:,j,i) = squeeze(ROI_trial(:,j,i));
            plotdata_smoothY(:,j,i) = mind_preprocess(plotdataY(:,j,i),5,0);
            plotdata_normY = plotdata_smoothY(:,j,i)/(max(plotdata_smoothY(:,j,i)));
            if i==1
                plot([0:100], plotdata_normY+(1*j),'Color', [.9 .5 .1])
            else
                plot([0:100], plotdata_normY+(1*j),'Color', [0 .5 0])
            end
            ylim([0 length(numtrial)+1]);
            xticks([0 33.33333 66.66667 100]);
            xticklabels({'0', '100' , '200', '300'})
            axis square
        end
    end
end

set(gcf, 'Position',[-2529 -290 1000 1000])


%% Plot the heat maps
load('cmap_heatmap3.mat');

figure; 
for k=1:length(randDoublet)
    curDoublet = randDoublet(k);
    ax{k} = subplot(5,5,k);
    roi1 = outputDoublets.saveAll_doublets(curDoublet).cells;
    roi1_all(k,:) = roi1;
    imagesc([mat2gray(out_Y.unfilteredAverageMaps(roi1(1)).averageMap) -1*mat2gray(out_Y.unfilteredAverageMaps(roi1(2)).averageMap)]')
    %imagesc([mat2gray(out_Y.pixelwise(roi1(1)).avgRealMap) -1*mat2gray(out_Y.pixelwise(roi1(2)).avgRealMap)]')
    colormap(cmap_heatmap3);
    xticks([.5 10.5 20.5 30.5]);
    xticklabels({'0', '100' , '200', '300'})
    caxis([-1 1]);
    set(gca,'box','off')
    axis square
end

set(gcf, 'Position',[-1529 -290 1000 1000])

