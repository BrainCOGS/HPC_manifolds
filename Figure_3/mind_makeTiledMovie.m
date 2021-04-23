function mind_makeTiledMovie(manifold3d, ROIactivities_thres, moviename)
% Need to run mind_scriptTiledFields to get the manifold3d and
% ROIactivities_thres variables

load('fieldwidth_colormap2.mat')

figure('Renderer', 'opengl', 'Position', [10 10 800 700]), clf;
for i=1:size(ROIactivities_thres,2)
    pos = 1+mod(i-1,5) + 6*(i-1-mod(i-1,5))/5;
    ax{i} = subplot(5,6,pos);
    ROInonzero = find(ROIactivities_thres(:,i));
    ROIzero    = boolean(ones(length(ROIactivities_thres),1));
    ROIzero(ROInonzero) = 0;

    hold on;
    h = scatter3(manifold3d(ROInonzero,1), manifold3d(ROInonzero,2), manifold3d(ROInonzero,3), 20, ROIactivities_thres(ROInonzero,i),'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8]);
    h2 = scatter3(manifold3d(ROIzero,1), manifold3d(ROIzero,2), manifold3d(ROIzero,3), 20, 'filled', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', .05);
    colormap(cmap_firingfield);
    caxis([0 3]);
    xlim([-0.07,0.07])
    ylim([-0.07,0.07])
    zlim([-0.07,0.07])

end

Link = linkprop([ax{:}],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});

for i=1:25
    set(gcf,'CurrentAxes',ax{i}); 
    set(ax{i},'CameraViewAngleMode','Manual')
    axis square
    axis off
    view(-45,45)
end


h = subplot(5,6,[6,12,18,24,30]);
[x,y] = meshgrid(0:0.4:1,0:0.1:5);
imagesc(flipud(y));
daspect([1 1 1])
set(h, 'position', [0.5 0.1 0.8 0.7] );
set(gca, 'XTick', [], 'XTickLabel', []) 
yticks([.5,26,51.5])
yticklabels({'≥2','1','0'})
a = get(gca,'yTickLabel');
set(gca,'yTickLabel',a,'fontsize',18)
text(0,-1.3,'\DeltaF/F','FontSize',20)
text(-58,-12,'Example manifold firing fields','FontSize',25)
ax1=gca;
set(ax1,'YAxisLocation','right');
set(ax1,'linewidth',3)
ax1.XAxis.Visible = 'off';
box off;

%%
moviename2 = strcat('./', moviename);
%video1 = VideoWriter('./test3','MPEG-4');
video1 = VideoWriter(moviename2,'MPEG-4');

open(video1);

ax{i} = subplot(5,6,1);
for i = 1:4:360 %360
   camorbit(4,0,'data',[0 1 0])
   drawnow
   writeVideo(video1,getframe(gcf));
end
for i = 1:4:360
   camorbit(4,0,'data',[1 0 0])
   drawnow
   writeVideo(video1,getframe(gcf));
end
for i = 1:4:360
   camorbit(4,0,'data',[0 0 1])
   drawnow
   writeVideo(video1,getframe(gcf));
end
close(video1);