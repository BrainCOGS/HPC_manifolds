function plot_manifold_figure(y, ontop, az, el, caxis_min, caxis_max, savepath, remove_axis)
% DESCRIPTION:
% This function saves a plot of manifold "y" with colorcoded variables "ontop", as
% viewed from azimuth "Az" and elevation "el" with the limits of the
% colorbar "caxis_min" and "caxis_max" to the location "savepath".
% remove_axis is a flag on whether axes are to be included in the pixelated graphics. 
% This is is for alignment-health-check. Usually false.
%
% OUTPUT:
% The saved file "myfigure_vec.pdf" contains axes and lines as vectors, but
% the scatterplot as pixel graphic. This makes for considerably smaller
% figures.
%
% EXAMPLE USE:
% plot_manifold_figure( randn(10000,3), randn(10000,1), 45, 35, 0, 10, '/Users/ms81/', false)


    close all;
    
    fig3flag = false;
    %plot_manifold_figure( all_manifold, smooth_variable, -63, 36, -9.85, 7.55, '/Users/ms81/', true);
    
    figure(1)
    set(gcf, 'units', 'points');
    figurePosition = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'points', 'PaperSize', [figurePosition(3) figurePosition(4)])
    set(gcf, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 figurePosition(3) figurePosition(4)]);

    rawAxis = axes(gcf, 'color', 'none', 'box', 'off', 'units', 'points');
    if fig3flag
        scatter3(rawAxis, y(:,1), y(:,2), y(:,3), 80, ontop, 'filled', 'MarkerFaceAlpha', .5);
        colorbar
        grid on;
        xlim([-0.07,0.07])
        ylim([-0.07,0.07])
        zlim([-0.07,0.07])
        axis square
        rotate3d on
    else
        scatter3(rawAxis, y(:,1), y(:,2), y(:,3), [], ontop, '.')
%         scatter3(rawAxis, y(:,1), y(:,2), y(:,3), [], ontop, '.')
        xlim([-0.14,0.14])
        ylim([-0.14,0.14])
        zlim([-0.14,0.14])
    end
    if remove_axis
        axis off;
    end
    daspect([1 1 1])
    caxis([caxis_min, caxis_max])
    view(az, el)
    set(colorbar,'visible','off')
    xl = xlim; yl = ylim; zl = zlim;
    print(gcf, [savepath, 'myfigure_raster.png'], '-dpng', '-r300')


    figure(3)
    set(gcf, 'units', 'points');
    set(gcf, 'PaperUnits', 'points', 'PaperSize', [figurePosition(3) figurePosition(4)])
    set(gcf, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 figurePosition(3) figurePosition(4)]);

    rasterAxis = axes(gcf, 'color', 'none', 'box', 'off', 'units', 'points');
    scatterAxis = axes(gcf, 'color', 'none', 'box', 'off', 'units', 'points');
    set(rasterAxis, 'position', [0 0 figurePosition(3) figurePosition(4)]);
    [A, ~, alpha] = imread([savepath, 'myfigure_raster.png']);
    if isempty(alpha)==1
        imagesc(rasterAxis, A);
    else
        imagesc(rasterAxis, A, 'alphaData', alpha);
    end
    axis(rasterAxis, 'off');
    hold on;
    s = scatter3(scatterAxis, 0, 0, 0, [], 0, 'w.')
    s.Marker = 'none';
    grid on;
    daspect([1 1 1])
    caxis([caxis_min, caxis_max])
    if fig3flag
        zticks([-0.05, 0, 0.05]);
    end
    view(az, el)
    set(colorbar,'visible','on')
    axis on;
    xlim(xl); ylim(yl); zlim(zl)
    print(gcf, [savepath, 'myfigure_vec.pdf'], '-dpdf')
    
end