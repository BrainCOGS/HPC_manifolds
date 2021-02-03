classdef FormatDefaults
   
  properties (Constant)

    %% Colors 
%     myBlue                = [57 133 255]./255;
    myBlue                = [0 123 255]./255;
    RchoiceCl             = FormatDefaults.myBlue;
    LchoiceCl             = [1 0 0];
    choiceCl              = [FormatDefaults.LchoiceCl; FormatDefaults.RchoiceCl];
%     choiceClMedium        = [227 102 102; 90 150 214; 224 188 224]/255;
    choiceClMedium        = [227 102 102; 90 150 214; 150 90 150]/255;
    choiceClLight         = bsxfun(@plus, 0.4*FormatDefaults.choiceCl, 0.6*[1 1 1]);
    choicePalette         = {@hot4, @cold4};        % L,R
    absDeltaCMap          = @purple2;
    perfBoutThCMap        = @green;
    correctCl             = [.75 .75 .75];
    errorCl               = [192 71 196]/255;
    errorClLight          = 0.5*FormatDefaults.errorCl + 0.5*[1 1 1];
%     layerCl               = [255 208 0; 158 85 36]/255;
%     layerClDark           = [196 160 0; 135 54 0]/255;
%     layerCl               = [0.45 0.45 0.45; 0.8 0.8 0.8];
    layerCl               = [204 191 224; 110 110 110]/256;
    layerClDark           = [154 110 219; 80 80 80]/256;
    lateralCl             = [139 196 164; 76 140 224; 63 209 217]/255;
%     areaCl                = [0 0 0; 120 11 11; 190 50 0; 204 100 0; 204 150 0; 181 181 0]/255;
%     areaCl                = [0 0 0; 204 126 204; 144 93 212; 53 115 232; 55 191 175; 130 194 47]/255;
    areaCl                = [0 0 0; 135 0 88; 190 0 0; 204 100 0; 204 150 0; 152 176 0]/255;
    annotationCl          = repmat([0.4; 0.7; 0.55], 1, 3);
    
    panelGray             = [.9 .9 .9];
    lightestGray          = [.95 .95 .95];
    lighterGray           = [.85 .85 .85];
    lightGray             = [.75 .75 .75];
    mediumGray            = [.6 .6 .6];
    darkerGray            = [.4 .4 .4];
    darkGray              = [.3 .3 .3];
    lightGreen            = [109 222 40]./255;
    mediumGreen           = [60 179 113]./255;
    darkGreen             = [0 .5 0];
    mediumLaurel          = [184 207 122]/255;
    darkLaurel            = [154 171 103]/255;
    lightRed              = [227 134 134]/255;
    brightRed             = [1 0 0];
    mediumRed             = [.7 0 0];
    darkRed               = [.5 0 0];
    lightBlue             = [151 193 252]./255;
    lightPurple           = [220 150 255]./255;
    darkBlue              = [0 65 179]/255;
    mediumPurple          = [200 112 230]./255;
    darkPurple            = [188 19  254]./255;
    mediumLavender        = [184 139 196]/255;
    darkLavender          = [156 118 166]/255;
    lightYellow           = [255 220 150]./255;
    darkYellow            = [252 194 0]./255;
    mediumOrange          = [235 153 0]/255;
    darkOrange            = [235 117 0]/255;
    isSigCl               = [.3 .3 .3];
    isNotSigCl            = [.85 .85 .85];
    animalEgCl            = {[195 102 13]./255,[35 113 20]./255,[0 0 0]};
    animalEgClLight       = {[193 128 33]./255,[50 125 30]./255,[.3 .3 .3]};
    axisCl                = [0 0 0];
    scatterCl             = [.3 .3 .3];
    genotypeDarkCl        = {[0 0 1],   [.3 .3 .3], [0 .5 0],   [.7 0 .5], [.96 .57 .13]};
    genotypeLightCl       = {[.5 .5 1], [.8 .8 .8], [.5 1 .5], [1 .5 1], [1 .8 .4]};
    mazeColorsReboot      = {[.75 .75 .75],[.5 .5 .5],[.25 .25 .25],[0 0 0],[1 .7 .2],[.6 0 0],[.8 0 0],[.8 .3 .3],[1 0 0],[.9 .2 .9],[.5 .2 .9],[.2 .5 1],[.2 .2 1],[0 0 1]};
    mazeColorsCondensed   = {[.75 .75 .75],[.5 .5 .5],[.25 .25 .25],[0 0 0],[.6 0 0],[1 .7 .2],[1 0 0],[.9 .2 .9],[.5 .2 .9],[0 0 .4],[.6 .6 1]};
    
    %% Lines, markers, and fonts
%     linewidthRegular      = .5;
%     linewidthThin         = .25;
%     linewidthThick        = .75;
    linewidthRegular      = 1;
    linewidthVeryThin     = 0.25;
    linewidthThin         = 0.5;
    linewidthThick        = 1.5;
    linewidthVeryThick    = 2;
    linewidthScaleBar     = 3;
    markerSize            = 5;
    markerSizeSmall       = 3;
    markerSizeLarge       = 7;
    markerSizeDot         = 10;
    titleFontSize         = 16;
    modelFontSize         = 20;
    fontSize              = 15;
    legendFontSize        = 12;
    barWidth              = 0.96;
    linestyleRef          = '-.';
    linestyleOverlay      = {'-', '--'};
    
    %% Annotations
    arrowSize             = 6;
    
  end
  
  
  methods (Static)

    % Interpolated colomap, using code from the sc package by Oliver Woodford
    % https://www.mathworks.com/matlabcentral/fileexchange/16233-sc-powerful-image-rendering
    function cmap = evidenceColors(brightOrDark, nColors)
      if nargin < 2
        nColors   = 512;
      end
      if brightOrDark
        choiceCl  = FormatDefaults.choiceClLight;
        bkgCl     = [1 1 1];
      else
        choiceCl  = FormatDefaults.choiceCl;
        bkgCl     = [1 1 1]*0;
      end
      cmap        = colormap_helper([ choiceCl(1,:)                         ...
                                    ; 0.3*mean(choiceCl,1) + 0.7*bkgCl      ...
                                    ; choiceCl(end,:)                       ...
                                    ], nColors);
%   colormap_helper(FormatDefaults.choiceClLight, 512);
%   colormap_helper([   0   0 100 1 ;   0   0 255 1 ;   0 150 255 1 ;   0 213 255 1 ; 200 255   0 1   ...
%                  ; 255 255   0 1 ; 255 150   0 1 ; 255   0   0 1 ; 100   0   0 1                   ...
%                  ]/256, 512);
%   colormap_helper([0 0 0 1; 80 80 80 1; 125 118 96 1; 171 145 109 1; 201 174 85 1; 224 188 70 1; 237 213 33 1; 255 234 0 1; 255 240 100 1; 245 245 130 1]/255, 512);
%   colormap_helper([100 100 100 1; 125 118 96 1; 161 149 112 1; 201 174 85 1; 224 188 70 1; 255 201 33 1; 255 234 0 1]/255, 512);
    end
    
    % Interpolated colomap, using code from the sc package by Oliver Woodford
    % https://www.mathworks.com/matlabcentral/fileexchange/16233-sc-powerful-image-rendering
    function cmap = strengthColors(nColors, minColor, maxColor)
      if nargin < 1
        nColors   = 512;
      end
      if nargin < 2
        minColor  = [0 0 0];
      end
      if nargin < 3
        maxColor  = FormatDefaults.darkOrange;
      end
      cmap        = colormap_helper([minColor; maxColor], nColors);
    end
    
    % Latex formatting for numbers with confidence intervals
    function str = numberWithCI(number, interval, scaleFactor)
      if ~exist('scaleFactor', 'var') || isempty(scaleFactor)
        scaleFactor = 1;
      end
      str           = sprintf ( '%.3g_{-%.2g}^{+%.2g}'                    ...
                              , scaleFactor *   number                    ...
                              , scaleFactor * ( number - interval(1) )    ...
                              , scaleFactor * ( interval(2) - number )    ...
                              );
    end
    
  end


end

