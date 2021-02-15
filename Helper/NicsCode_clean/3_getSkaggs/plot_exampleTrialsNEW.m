function out_examples = plot_exampleTrialsNEW(fname, skaggs_out, whichTrials)

% taskType is either 'towers' or 'alternation'
% trialType is 'goodTrials' or 'keepTrials'
% shuffleName is of the format 'E66_Y_[11_4]_%d.mat'
% in general, whichTrials should be set to be exactly the same as the
% argins from skaggs_out, i.e. skaggs_out.argins.trialsAverageMap
% but can also set to something else

% load(fname)
% if strcmp(trialType, 'keepTrials')
%     keepTrials = find([score.trial.mainTrial]==1 & [score.trial.goodQuality]==1);
% elseif strcmp(trialType, 'goodTrials')
%     keepTrials = find([score.trial.goodQuality]==1);
% end
switch skaggs_out.argins.taskType
    case 'towers'
        config.trialLength = 330;
    case 'alternation'
        config.trialLength = 370;
end

out.config = config;

%%

numROI = length(skaggs_out.ROIs);

nic_output_POS = extractVariables('all', 4, whichTrials, 2, config.trialLength, 0, 5, skaggs_out.argins.DFFthresholding, fname,'none',skaggs_out.argins.taskType,1,1);
behavioralVariables = nic_output_POS.behavioralVariables;
trialn   = behavioralVariables.Trial;
trials = unique(trialn);
numtrials = length(trials);
%numtrials = length(skaggs_out.trials);

pfInfo_this = skaggs_out.skaggsMetric.skaggs_real;

[~, maxint] = sort(pfInfo_this,'descend');

DFF_trials = reshape(nic_output_POS.ROIactivities,[config.trialLength+1,numtrials,numROI]);



layout       = [ 1 2 3 4 5 6 7 8 9 10];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;
% figure
% numplot = 6;
% c=1;
for i=1:length(layout)
    % subplot(1,numplot,c);
    iPanel      = iPanel + 1;
    axs         = fig.panel(iPanel);
    DFF_trials1 = DFF_trials(:,:,maxint(i))';
    %     for i=1:size(DFF_trials1,1)
    %         DFF_norm(i,:) = mat2gray(DFF_trials1(i,:));
    %     end
    
    % imagesc(DFF_trials(:,:,maxint(i))');
    
    % imagesc(DFF_norm);
    imagesc(mat2gray(DFF_trials1));
    
    xticks([30.5 130.5 230.5 330.5]);
    xticklabels({'0','100','200','300'});
    xlabel('Position (cm)');
    ylabel('Trial #');
    set(gca,'box','off')
    set(gca, 'FontName', 'Arial')
    
    % c=c+1;
    
end

% set(gcf, 'Position', [-1200, -550, 800, 1200])

out.DFF_trials = DFF_trials;
