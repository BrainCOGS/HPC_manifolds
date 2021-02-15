function outputMindPlotter = mind_plotManifoldGradients(outMind, fname, taskType, togglePlot)

choiceORviewangle = 'viewangle';

% togglePlot = 1 means to plot it
dimEmbed = 3;

sample = outMind.dataDFF;
pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
manifold3d = outMind.dat.allembed(outMind.config_input.mindparameters.embed.d==dimEmbed).f2m.map.transform(pca_coords);

if strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [11 4], fname,'none','towers',1,1);
elseif strcmp(taskType, 'alternation')==1 || strcmp(taskType, 'Alternation')==1
    nic_output = extractVariables('all', 2, 'goodTrials', 2, 0, 0, 5, [11 4], fname,'none','alternation',1,1);
end
behavioralVariables = nic_output.behavioralVariables;
ROIactivities = nic_output.ROIactivities;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part if correctOnly
% choiceCorrect = nic_output.behavioralVariables.ChoiceCorrect;
% ROIactivities = nic_output.ROIactivities;
% ROIactivities = ROIactivities(choiceCorrect==1,:);
% behavioralVariables = behavioralVariables(choiceCorrect==1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if outMind.argins.downSample==0
%     downSample = 1;
% else
%     downSample = outMind.argins.downSample;
% end
%behavioralVariables = behavioralVariables(1:outMind.argins.downSample:end,:);
behavioralVariables = behavioralVariables(outMind.Datarange,:);

if strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
    smooth_E = mind_smoothDimensions(manifold3d,behavioralVariables.Evidence,20);
end
smooth_Y = mind_smoothDimensions(manifold3d,behavioralVariables.Position,20);
smooth_V = mind_smoothDimensions(manifold3d,rad2deg(behavioralVariables.ViewAngle),20);
Choice1  = behavioralVariables.Choice;

%ROIactivities = ROIactivities(1:outMind.argins.downSample:end,:);
ROIactivities = ROIactivities(outMind.Datarange,:);
meanFR = mean(ROIactivities, 2);





outputMindPlotter.manifold3d = manifold3d;
outputMindPlotter.dimEmbed = dimEmbed;
outputMindPlotter.smooth_Y = smooth_Y;
if strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
    outputMindPlotter.smooth_E = smooth_E;
end
outputMindPlotter.Choice1 = Choice1;

if togglePlot==1
    figure;
    ax(1) = subplot(1,3,1);
    hold on;
    scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3), 20, smooth_Y,'filled','MarkerFaceAlpha',.5);
    xlabel('Dim 1')
    ylabel('Dim 2')
    zlabel('Dim 3')
    title('Position (smoothed)')
    %axis tight
    xlim([-0.07,0.07])
    ylim([-0.07,0.07])
    zlim([-0.07,0.07])
    axis square
    colorbar
    grid on;
    
    
    if strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
        ax(2) = subplot(1,3,2);
        hold on;
        scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3), 20, smooth_E,'filled','MarkerFaceAlpha',.5);
        xlabel('Dim 1')
        ylabel('Dim 2')
        zlabel('Dim 3')
        title('Evidence (smoothed)')
        %axis tight
        xlim([-0.07,0.07])
        ylim([-0.07,0.07])
        zlim([-0.07,0.07])
        axis square
        colorbar
        grid on;
        
    end
    
    
    ax(3) = subplot(1,3,3);
    hold on;
    if strcmp(choiceORviewangle, 'choice')
        scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3), 20, Choice1,'filled','MarkerFaceAlpha',.5);
        title('Choice')
    elseif strcmp(choiceORviewangle, 'viewangle')
        scatter3(manifold3d(:,1), manifold3d(:,2), manifold3d(:,3), 20, smooth_V,'filled','MarkerFaceAlpha',.5);
        title('View Angle')
        caxis([-30, 30]);
    end
    xlabel('Dim 1')
    ylabel('Dim 2')
    zlabel('Dim 3')
    
    colorbar
    grid on;
    %axis tight
    xlim([-0.07,0.07])
    ylim([-0.07,0.07])
    zlim([-0.07,0.07])
    axis square
    
    
    
    rotate3d on
    
    set(gcf, 'Position', [100, 340, 1300, 420])
    
    hlink = linkprop([ax],{'CameraPosition','CameraUpVector'});
    
    % For E65
    %view(-18, 22);
    %view(-142, 21);
    % for optimized
    %view(-51, -35);
    % for ml500_lmf1
    view(-63,36)
    
    % For E47
    %view(11, -12.5);
    outputMindPlotter.hlink = hlink;
end
