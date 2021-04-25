function plot_getSkaggs_summarySLIM(getSkaggsOutput, selectedROIs, numberPlot, dimPlot)

% Creates figures from the output of the Mutual Information analysis.
%
% Sample call:
% plot_getSkaggs_summary(getSkaggsOutput, getSkaggsOutput.skaggsMetric.sigROIs, 25, [5, 5])
%
% ***** INPUTS *****
% getSkaggsOutput - output structure from getSkaggs()
% selectedROIs - vector of ROI labels to plot
% numberPlot - double specifying how many plots to make. should be equal to
%               or less than dimPlot(1)*dimPlot(2)
% dimPlot - two-length vector of the subplot dimensions. ex: [5, 5];
%
% *** OUTPUTS: ***
% figHandles - a cell array containing figure handles to the following
% figures:
%
% figure 1 - Histogram of MI ROIs in real data vs shuffle tests
% figure 2 -
% figure 3 -
% figure 4 -
% figure 5 -
% figure 6 -

%% Prepare data
% save relevant fields of MI analysis output. here, the fields from
% pixelwise will be used instead of unfilteredAverageMaps because
% visualization  will use the data that has been smoothed (if applied) with
% a Gaussian filter.

pixelwise = getSkaggsOutput.pixelwise;
pX = getSkaggsOutput.pX;
skaggsMetric = getSkaggsOutput.skaggsMetric;

% check that the number of dimensions is equal to 2
if length(getSkaggsOutput.argins.dimensions) ~= 2
    error('Error: This function only works for 2D analyses.');
end


%% Plot each field using surf for real lambda(x)
% POSITION AND EVIDENCE HAVE BEEN HARD-CODED HERE
figure;
for i = 1:numberPlot

    ind = selectedROIs(i); % get ROI index for pixelwise
    
    subplot(dimPlot(1), dimPlot(2), i); % set subplot index
    
    z = pixelwise(ind).avgRealMap; % real lambda(x)
    z(z < 0) = 0; % for visualization, set negative lambda(x) to zero
    
    surf(z, 'faceColor', 'interp'); % create surface using interpolation
    set(gca,'color','none'); % set color to none
    
    
    hold on; 

    ylim([0 size(pX, 1)]);
    xlim([0 size(pX, 2)]);
    zlim([0 max(max(z))+.1]);
    ylabel(sprintf('ROI %d', ind));
    yticks([6 16 26]); 
    xticks([0 10 20 30]); 
    
    if i==1
        xticklabels({'0' , '100' , '200', '300'})
        yticklabels({'-10' , '0', '10'})
    else
        xticklabels({});
        yticklabels({});
    end
    view(-40,70);
    
    % axis square
    if i == 2
        title('Selected ROIs Average Map');
    end
end
clear i
set(gcf, 'Position',[-2268 -340 1300 1000])



%% avg signif map, but each ROI is its own color

colors = distinguishable_colors(numberPlot);
% colors = {[0.3 0.6 0.3], [0.6 0.3 0.3], [0.3 0.3 0.6], [0.6 0.3 0.6]};
figure;
% colorCounter = 4;
for i = 1:numberPlot
    
    hSurface = surf(pixelwise(selectedROIs(i)).realSignifMap);
    set(hSurface, 'FaceColor', colors(i,:), 'faceAlpha', 0.2);
    % set(hSurface,'LineWidth',.1)
    hold on
    xticks([0 10 20 30])
    xticklabels({'0', '100', '200', '300'});
    xlabel('Position (cm)');
    yticks([6 16 26])
    yticklabels({'-10','0','10'});
    ylabel('Evidence');
    % drawnow
    %     colorCounter = colorCounter + 1;
    %     if colorCounter == size(colors, 2) + 1
    %         colorCounter = 1;
    %     end
    
end
title('Individual ROI Fields')
set(gca,'color','none')
view(-45, 60);


%% Plot just the single example
figure;
i=1;
ind = selectedROIs(i); % get ROI index for pixelwise

%subplot(dimPlot(1), dimPlot(2), i); % set subplot index

z = pixelwise(ind).avgRealMap; % real lambda(x)
z(z < 0) = 0; % for visualization, set negative lambda(x) to zero

surf(z, 'faceColor', 'interp'); % create surface using interpolation
set(gca,'color','none'); % set color to none


hold on;

ylim([0 size(pX, 1)]);
xlim([0 size(pX, 2)]);
zlim([0 max(max(z))+.1]);
ylabel(sprintf('ROI %d', ind));
yticks([6 16 26]);
xticks([0 10 20 30]);


xticklabels({'0' , '100' , '200', '300'})
yticklabels({'-10' , '0', '10'})

view(-40,70);

% axis square
if i == 2
    title('Selected ROIs Average Map vs Shuffle Shroud (surface)');
end
view(-45, 60);


