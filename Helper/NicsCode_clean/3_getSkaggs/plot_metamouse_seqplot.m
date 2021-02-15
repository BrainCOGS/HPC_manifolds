function [metaTrain, metaTest] = plot_metamouse_seqplot(metamouse)
% 
% Description: plots 1-D metamouse sequence plots for the metamouse. This
% function uses the output structure from generateMetamouse, which may use
% both left- and right-choice trials, or split the trials into left- and
% right-choice only. The goal is to concatenate the sequence plots across
% animals and resort their sequence plots for the metamouse.
%
% Sample Call: [metamouseTraining, metamouseTesting] = plot_metamouse_seqplot(metamouse_output)
% 
% ***** INPUTS *****
% metamouse - output structure from generateMetamouse()
% 
% ***** OUTPUTS *****
% metaTrain - matrix of sorted ROIs' mean DFF in training set across animals
% metaTest - matrix of sorted ROIs' mean DFF in testing set across animals
% a figure showing the individual animals' sequence plots and metamouse

% preallocate output variables
metaTrain = []; % training set
metaTest = [];  % testing set
metaAll = [];   % no CV set

% the metamouse strcture did *not* split between left- and right-choice
% trials
if all(isfield(metamouse, {'LC', 'RC'})) == 0
    % concatenate the individual animals' sequence plots into a metamouse
    % sequence plot
    for i = 1:length(metamouse)
        metaTrain = [metaTrain; ...
            metamouse(i).getSkaggs.sequencePlot.CV.matrixTrain]; %#ok<*AGROW>
        metaTest = [metaTest; ...
            metamouse(i).getSkaggs.sequencePlot.CV.matrixTest];
        metaAll  = [metaAll; ...
            metamouse(i).getSkaggs.sequencePlot.noCV.sortedMatrix];
    end
    
    % for the metamouse training set, find the new sorting sequence based
    % on location of the max mean DFF value for each ROI
    [~, metaMaxTrain] = max(metaTrain, [], 2);
    [~, metaSortTrain] = sort(metaMaxTrain);
    
    % repeat this process for the no CV set
    [~, metaMaxAll] = max(metaAll, [], 2);
    [~, metaSortAll] = sort(metaMaxAll);
    
    %% Plot individual animals' plots in one figure
    % create a new figure to plot the individual animals' sequence plots
    % for the training, testing, and no CV sets
    figure;
    figDim2 = length(metamouse)+1;
    for i = 1:length(metamouse) % loop through animals
        % plot training set
        subplot(3,figDim2, i);
        plotPanels(metamouse(i).getSkaggs.sequencePlot.CV.matrixTrain, ...
            metamouse(i).getSkaggs, [metamouse(i).animalID, ' Train Set'], 'Neuron #');
        
        % plot testing set
        subplot(3,figDim2, i+figDim2);
        plotPanels(metamouse(i).getSkaggs.sequencePlot.CV.matrixTest, ...
            metamouse(i).getSkaggs, [metamouse(i).animalID, ' Test Set'], 'Neuron #');
        
        % plot no CV set
        subplot(3,figDim2, i+(figDim2*2));
        plotPanels(metamouse(i).getSkaggs.sequencePlot.noCV.sortedMatrix, ...
            metamouse(i).getSkaggs, [metamouse(i).animalID, ' All Trials Set'], 'Neuron #');
        
    end
    
    %% Plot the metamouse sequence plot
    
    % resort the metamouse training and testing sets based on the training
    % sort index
    metaTrain = metaTrain(metaSortTrain, :);
    metaTest = metaTest(metaSortTrain, :);
    % resort the no CV set based on the no CV sorting index
    metaAll = metaAll(metaSortAll, :);
    
    % add the metamouse plots to the figure
    % plot the training set
    subplot(3,figDim2,figDim2);
    plotPanels(metaTrain, metamouse(1).getSkaggs, 'Meta Train', 'Neuron #');
    
    % plot the testing set
    subplot(3,figDim2,2*figDim2);
    plotPanels(metaTest, metamouse(1).getSkaggs, 'Meta Test', 'Neuron #');
    
    % plot the no CV set
    subplot(3,figDim2,3*figDim2);
    plotPanels(metaAll, metamouse(1).getSkaggs, 'Meta All', 'Neuron #');
    
    % adjust figure position
    set(gcf, 'Position', [-2500, 50, 1900, 700])
    
else
    % in this case, the metamouse analysis *was* split between left- and
    % right-choice trials. an ROI that was significant in *either* the 
    % left-choice (LC) or right-choice (RC) analysis 'preferred' that 
    % choice, whereas an ROI that was significant in *both* LC and RC
    % analyses was 'non-specific'. Look at the sequence plots for the
    % preferred ROIs on their preferred and non-preferred trials
    
    RC = [metamouse.RC];            % right-choice trial analysis
    LC = [metamouse.LC];            % left-choice trial analysis
    numAnimals = length(metamouse); % number of animals
    
    % preallocate matrices
    mRC_pref = [];                  % RC ROIs on preferred (RC) trials
    mRC_nonPref = [];               % RC ROIs on non-preferred (LC) trials
    mLC_pref = [];                  % LC ROIs on preferred (LC) trials
    mLC_nonPref = [];               % LC ROIs on non-preferred (RC) trials
    mNonSpecific_LC = [];           % non-specific ROIs on LC trials
    mNonSpecific_RC = [];           % non-specific ROIs on RC trials
    
    %% from each animal, extract the relevant matrices
    for i = 1:numAnimals
        % RC_ROIs = RC(i).skaggsMetric.sigROIs;
        % LC_ROIs = LC(i).skaggsMetric.sigROIs;
        RC_ROIs = RC(i).sequencePlot.ROIsPassedCV; % RC ROI indices
        LC_ROIs = LC(i).sequencePlot.ROIsPassedCV; % LC ROI indices
        % non-specific ROIs are common to both RC and LC ROIs
        nonSpecific_ROIs = unique([RC_ROIs(ismember(RC_ROIs, LC_ROIs)), ...
            LC_ROIs(ismember(LC_ROIs, RC_ROIs))]);
        % delete non-specific ROIs from the LC and RC ROIs lists
        RC_ROIs = RC_ROIs(~ismember(RC_ROIs, nonSpecific_ROIs));
        LC_ROIs = LC_ROIs(~ismember(LC_ROIs, nonSpecific_ROIs));
        
        % get the data for the RC preferring ROIs on preferred (RC) and
        % non-preferred (LC) trials
        mRC_pref    = [mRC_pref, ...
            RC(i).unfilteredAverageMaps(RC_ROIs).averageMap]; %#ok<*AGROW>
        mRC_nonPref = [mRC_nonPref, ...
            LC(i).unfilteredAverageMaps(RC_ROIs).averageMap];
        
        % follow a similar procedure for LC preferring ROIs
        mLC_pref    = [mLC_pref, ...
            LC(i).unfilteredAverageMaps(LC_ROIs).averageMap];
        mLC_nonPref = [mLC_nonPref, ...
            RC(i).unfilteredAverageMaps(LC_ROIs).averageMap];
        
        % for non-specific ROIs, get their data for LC and RC trials
        % separately
        mNonSpecific_LC = [mNonSpecific_LC, ...
            LC(i).unfilteredAverageMaps(nonSpecific_ROIs).averageMap];
        mNonSpecific_RC = [mNonSpecific_RC, ...
            RC(i).unfilteredAverageMaps(nonSpecific_ROIs).averageMap];
    end
    
    %% Normalize each ROI by row. 
    % because the number of ROIs in each list (RC, LC, and non-specific 
    % ROIs), there's a different loop for each list. 
    
    % RC ROIs
    for i = 1:size(mRC_pref, 2)
        % preferred trials
        mRC_pref(~isnan(mRC_pref(:,i)), i) = ...
            mat2gray(mRC_pref(~isnan(mRC_pref(:,i)),i));
        % non-preferred trials
        mRC_nonPref(~isnan(mRC_nonPref(:,i)),i) =  ...
            mat2gray(mRC_nonPref(~isnan(mRC_nonPref(:,i)),i));
    end
    
    % LC ROIs
    for i = 1:size(mLC_pref, 2)
        % preferred trials
        mLC_pref(~isnan(mLC_pref(:,i)), i) = ...
            mat2gray(mLC_pref(~isnan(mLC_pref(:,i)),i));
        % non-preferred trials
        mLC_nonPref(~isnan(mLC_nonPref(:,i)),i) =  ...
            mat2gray(mLC_nonPref(~isnan(mLC_nonPref(:,i)),i));
    end
    
    % non-specific ROIs
    for i = 1:size(mNonSpecific_LC, 2)
        % LC trials
        mNonSpecific_LC(~isnan(mNonSpecific_LC(:,i)),i) = ...
            mat2gray(mNonSpecific_LC(~isnan(mNonSpecific_LC(:,i)),i));
        % RC trials
        mNonSpecific_RC(~isnan(mNonSpecific_RC(:,i)),i) = ...
            mat2gray(mNonSpecific_RC(~isnan(mNonSpecific_RC(:,i)),i));
    end
    
    %% Sort the matrices
    % for the LC and RC preferring ROIs, sort the ROIs in the preferred and
    % non-preferred trial sets based on the ROIs' maximum mean DFF on the
    % preferred trial set
    
    % LC ROIs
    [~, LC_maxInd_pref]     = max(mLC_pref); 
    [~, LC_sortInd_pref]    = sort(LC_maxInd_pref); % sort index
    
    mLC_pref    = mLC_pref(:, LC_sortInd_pref)'; 
    mLC_nonPref = mLC_nonPref(:, LC_sortInd_pref)'; 

    % follow a similar procedure for the RC ROIs
    [~, RC_maxInd_pref]     = max(mRC_pref);
    [~, RC_sortInd_pref]    = sort(RC_maxInd_pref); % sort index
    
    mRC_pref    = mRC_pref(:, RC_sortInd_pref)';
    mRC_nonPref = mRC_nonPref(:, RC_sortInd_pref)';
    
    % for non-specific ROIs, sort the LC and RC sets independently of each
    % other based on the location of the ROIs' maximum mean DFF value
    % LC trials only
    [~, NS_maxInd_LC] = max(mNonSpecific_LC);
    [~, NS_sortInd_LC] = sort(NS_maxInd_LC);
    mNonSpecific_LC = mNonSpecific_LC(:, NS_sortInd_LC)';
    
    % RC trials only
    [~, NS_maxInd_RC] = max(mNonSpecific_RC);
    [~, NS_sortInd_RC] = sort(NS_maxInd_RC);
    mNonSpecific_RC = mNonSpecific_RC(:, NS_sortInd_RC)';
    
    
    %% Plot the results
    % show the sequence plots for the LC- and RC-preferring ROIs on their
    % preferred and non-preferred trial sets. also, plot the non-specific
    % ROIs on the LC and RC trials
    
    figure; % create new figure
    % colormap(parula);
    
    % LC ROIs on LC trials
    subplot(3,2,1); 
    plotPanels(mLC_pref, metamouse(1).LC, ...
        'Left Choice (LC) Trials', 'LC-Specific ROIs');
    colorbar;
    
    % LC ROIs on RC trials
    subplot(3,2,2);
    plotPanels(mLC_nonPref, metamouse(1).LC, ...
        'Right Choice (RC) Trials', '');
    colorbar;
    
    % RC ROIs on LC trials
    subplot(3,2,3);
    plotPanels(mRC_nonPref, metamouse(1).RC, '', 'RC-Specific ROIs');
    colorbar;
    
    % RC ROIs on RC trials
    subplot(3,2,4);
    plotPanels(mRC_pref, metamouse(1).RC,'','');
    colorbar;
    
    % non-specific ROIs on LC trials
    subplot(3,2,5);
    plotPanels(mNonSpecific_LC, metamouse(1).LC, '', 'Choice Non-specific ROIs');
    colorbar;
    
    % non-specific ROIs on RC trials
    subplot(3,2,6);
    plotPanels(mNonSpecific_RC, metamouse(1).RC, '','');
    colorbar;
    
end
end

%% Function to plot individual panels of figure
function plotPanels(plotData, skaggsData, plotTitle, yText)
% Description: plots individual panels of function figure with appropriate
% title and label
%
% Sample Call: 
% plotPanels(mRC_pref, metamouse(1).RC,'','');
%
% ***** INPUTS *****
%
% plotData - matrix of data to plot with imagesc()
% skaggsData - analysis results from getSkaggs
% plotTitle - character vector for subplot title
% yText - character vector for subplot y-label


imagesc(plotData); % plot data
%xticks(''); yticks('');
title(plotTitle); % label title
xlabel(skaggsData.argins.dimensions(1)); % label x-label

% label the x-ticks if the dimension analyzed is Evidence or Position. if
% so, find the bin edges and label the xticks
switch skaggsData.argins.dimensions{1}
    case 'Evidence'
        binEdges = skaggsData.argins.binEdges{1};
        xticks([find(binEdges==-10) find(binEdges==0) find(binEdges==10)]);
        xticklabels({'-10', '0', '10'});
    case 'Position'
        binEdges = skaggsData.argins.binEdges{1};
        xticks([find(binEdges==0)-.5 find(binEdges==100)-.5 ...
            find(binEdges==200)-.5 find(binEdges==300)-.5]);
        xticklabels({'0', '100', '200', '300'});
end
% change figure settings
set(gca,'box','off')
set(gca, 'FontName', 'Arial')
% label y-axis
ylabel(yText);
end