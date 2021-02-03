% PANELEDFIGURE  Makes a new standardized figure with multiple panels for publications
% Runnable example:
%
%   fig = PaneledFigure([1 4 6; 2 4 6; 3 5 7],'small',[]);
%   xlabel(fig.panel(1),'Hello World'); 
%   ylabel(fig.panel(1),'Big plans');
%   fig.hint();
%
classdef PaneledFigure < handle
  
  %________________________________________________________________________
  %------- Constants
  properties (Constant)
    
    %% Controls output
    EXPORT_FORMAT     = {'.png', '.pdf'}              % always have png or some kind of Matlab readable image first for diffs
    EXPORT_OPTIONS    = {'-r300', '-painters'}        % default export_fig() options
    MONITOR           = PaneledFigure.largestMonitor()
    PANEL_LABEL       = double('A')
    MAX_DISP_ITEMS    = 100;
    MAX_DISP_LENGTH   = 20;                           % for any one dimension of a tensor
    NUMBER_PRECISION  = 6;
    
    %% Controls sizes of panels and margins, in pixels
    PIX_TO_INCHES     = 1/72                          % following Macs
    SIZE_SETTINGS     = {'tiny','smaller','small', 'medium', 'tall'}
    PANEL_SIZE        = [ 120 120; 160 160; 200 200; 300 300 ; 350 271 ]   % height, width
    FONTSIZE_FACTOR   = [ .75 1 1 1 1]
    FIGURE_MARGIN     = [ 10 10 10 10 ]               % top,bottom,left,right
    PANEL_SPACING     = [ 12 15 ]                     % vertical,horizontal
    AXES_MARGIN       = [ 10 40 40 10 ]               % top,bottom,left,right
    COPY_BORDER       = [ 0 0; 0 0 ]
    TICK_FONTSIZE     = 13
    LABEL_FONTSIZE    = 16
    PANEL_FONTSIZE    = 20
    
  end
  
  %------- Private data
  properties (GetAccess = protected)
    axs
    panelGrid
    gridIndex
  end
  
  %------- Public data
  properties (SetAccess = protected)
    size
    axesMargin

    fig
    groupLabel
  end
  properties
    firstLabel    = PaneledFigure.PANEL_LABEL
    newGroup
  end
  
  %________________________________________________________________________
  methods (Static)
    
    %----- This flag can be set to turn off time-costly functions i.e. code version tracking and export of individual panels
    function yes = developmentMode(newValue)
      
      persistent        inDevelopement;
      if isempty(inDevelopement)
        inDevelopement  = false;
      end
      
      %% If input is provided, set the current mode to the desired value, otherwise return the current mode
      if nargin > 0
        inDevelopement  = newValue;
      end
      yes               = inDevelopement;
      
    end
    
    %----- Set/retrieve global format settings; see comments within the code for how these are defined
    function value = formatOptions(what, name, value)
      
      %% The fields of the format struct defines display items that can be configured
      persistent              format;
      if isempty(format)
        format                = struct();
        format.panelLabel     = struct( 'FontWeight'  , 'bold'    ...
                                      );
      end
      
      if ~isfield(format, what)
        error('PaneledFigure:formatOptions', 'what must be one of the predefined configurables: %s', strjoin(fieldnames(format),', '));
      end
  
      %% Either set or retrieve the desired setting depending on whether the user provides a set value
      if nargin > 2
        %% If a value is provided, set it, creating a new entry if neccessary
        format.(what).(name)  = value;
      elseif nargin > 1
        %% If no value is provided, check that the setting exists and return it
        if ~isfield(format.(what), name)
          error('PaneledFigure:formatOptions', '''%s'' has not been configured for %s.', name, what);
        end
        value                 = format.(what).(name);
      else
        %% If no specific name is provided, return all settings associated to this display item type
        value                 = fieldnames(format.(what))';
        value(end+1,:)        = struct2cell(format.(what));
        value                 = value(:)';
      end
      
    end
    

    %----- Helper function to resize panels in a given direction (horizontal or vertical)
    function axs = distributePanels(axs, iPos, iSize, index, relativeSize)
      %% Also allow for uneven widths/heights based on the relativeSize input parameter
      if nargin < 5 || isempty(relativeSize)
        relativeSize        = ones(size(index));
      end
      relativeSize          = relativeSize(:)';
      
      %% Get current panel positions in ascending order
      panelPos              = accumfun(1, @(x) get(axs(x), 'Position'), index);
      [~,iOrder]            = sort(panelPos(:,iPos));
      index                 = index(iOrder);
      panelPos              = panelPos(iOrder,:);
      iOrder                = iOrder(1:numel(relativeSize));
      relativeSize          = relativeSize(iOrder - min(iOrder) + 1);
      
      %% Decide on an inter-panel offset that maintains the extent of the figure originally covered 
      [numCols, colPanel, colIndex]         ...
                            = uniquetol(panelPos(1:numel(relativeSize),iPos), 1e-4);
      numCols               = numel(numCols);
      totalExtent           = sum(panelPos(colPanel(end),[iPos,iSize])) - panelPos(1,iPos);
      if numCols > 1
        posOffset           = (totalExtent - sum(panelPos(colPanel,iSize))) / (numCols-1);
      else
        posOffset           = 0;
      end
      
      %% Redistribute sizes so that the summed size remains the same but with the desired relative sizes of panels
      totalSize             = sum(panelPos(end,[iPos,iSize])) - panelPos(1,iPos) - posOffset * (numCols-1);
      targetSize            = totalSize .* relativeSize / sum(relativeSize(colPanel));
      
      %% If the number of relativeSize targets is less than the number of panels, delete panels
      if numel(relativeSize) < numel(index)
        range               = index(numel(relativeSize)+1:end);
        delete(axs(range));
        axs(range)          = gobjects(1);
        panelPos(numel(relativeSize)+1:end,:) = [];
        index               = index(1:numel(relativeSize));
      end
      
      %% Recompute and apply panel positions
      panelPos(:,iSize)     = targetSize;
      for iPanel = 2:numel(index)
        panelPos(iPanel,iPos) = panelPos(1,iPos) + sum(targetSize(1:colIndex(iPanel)-1)) + posOffset*(colIndex(iPanel)-1);
      end
      for iPanel = 1:numel(index)
        set(axs(index(iPanel)), 'Position', panelPos(iPanel,:));
      end
    end

    %----- Helper function to compute an offset label location for either linear or logarithmically scaled axes
    function loc = offsetLocation(offset, range, index, inches, scale)
      if strcmpi(scale, 'log')
        range = log(range);
        loc   = exp( range(index) + offset*diff(range)/inches );
      else
        loc   = range(index) + offset*diff(range)/inches;
      end
    end
    
    %----- Helper function to convert a cell array to printable strings
    function [info, infoSize, colWidth, tableWidth] = formatCells(info)
      
      %% Get string representation of all elements
      info        = info(:,:,:);
      infoSize    = [size(info), 1];
      for iElem = 1:numel(info)
        if isnumeric(info{iElem})
          info{iElem} = sprintf('%.*g', PaneledFigure.NUMBER_PRECISION, info{iElem});
        end
      end
      elemWidth   = cellfun(@numel, info);
      colWidth    = permute(elemWidth, [2 1 3]);
      colWidth    = max(reshape(colWidth, size(colWidth,1), []), [], 2);
      tableWidth  = sum(colWidth) + 3*numel(colWidth);

    end
    
  end
  
  %________________________________________________________________________
  methods

    %----- Constructor
    % nPanels   : [nRows, nCols] --OR e.g.-- [ 1 3 4
    %                                        ; 2 3 4
    %                                        ; 5 6 7 ]
    function obj = PaneledFigure(nPanels, panelSize, figureMargin, panelSpacing, axesMargin, varargin)

      %% Default arguments
      if nargin < 3 || isempty(figureMargin)
        figureMargin        = PaneledFigure.FIGURE_MARGIN;
      end
      if nargin < 4 || isempty(panelSpacing)
        panelSpacing        = PaneledFigure.PANEL_SPACING;
      end
      if nargin < 5 || isempty(axesMargin)
        axesMargin          = PaneledFigure.AXES_MARGIN;
      end
      obj.axesMargin        = axesMargin;
      
      
      %% Parse input 
      if numel(nPanels) == 2 && isvector(nPanels)
        obj.size            = nPanels;
        panelGrid           = reshape(1:prod(nPanels), flip(nPanels))';
      else
        obj.size            = size(nPanels);
        panelGrid           = nPanels;
      end
      if ischar(panelSize)
        panelHW             = PaneledFigure.PANEL_SIZE(strcmpi(PaneledFigure.SIZE_SETTINGS, panelSize), :);
        fontMultiplier      = PaneledFigure.FONTSIZE_FACTOR(strcmpi(PaneledFigure.SIZE_SETTINGS, panelSize));
      else
        panelHW             = panelSize;
        fontMultiplier      = 1;
      end
      obj.panelGrid         = panelGrid;
      
      %% Compute figure and grid spacings (all units here in pixels)
      panelStride           = panelSpacing + panelHW;
      figureHeight          = sum(figureMargin([2 1])) + obj.size(1)*panelStride(1) - panelSpacing(1);
      figureWidth           = sum(figureMargin(3:4)) + obj.size(2)*panelStride(2) - panelSpacing(2);
      position              = [ 50, 50, figureWidth, figureHeight ];
      gridY                 = figureMargin(2) + (obj.size(1)-1:-1:0) * panelStride(1);
      gridX                 = figureMargin(3) + (0:obj.size(2)-1)    * panelStride(2);
      
      %% Create figure on the desired monitor (using inches here to retain Mac/Windows similarity)
      screenSize            = get(0, 'MonitorPosition');
      monitor               = PaneledFigure.MONITOR;
      if monitor < 0
        monitor             = size(screenSize, 1) + monitor+1;
      end
      screenSize            = screenSize(min(monitor,end), :);
      position(1:2)         = position(1:2) + screenSize(1:2);
      obj.fig               = figure( 'Units'             , 'inches'                                  ...
                                    , 'Position'          , position*PaneledFigure.PIX_TO_INCHES      ...
                                    , 'Resize'            , 'off'                                     ...
                                    , 'Color'             , [1 1 1]                                   ...
                                    , 'PaperPositionMode' , 'auto'                                    ...
                                    , varargin{:}                                                     ...
                                    );
      
      
      %% Create axes aligned to the specified grid
      nPanels               = max(panelGrid(:));
      obj.axs               = gobjects(1,nPanels);
      obj.gridIndex         = nan(1,nPanels);
      for iPanel = 1:nPanels
        if ~any(panelGrid(:) == iPanel)
          continue;
        end
        
        %% panelLoc(first/count, row/col)
        panelLoc            = {any(panelGrid == iPanel, 2), any(panelGrid == iPanel, 1)};
        panelSpan           = [ find(panelLoc{1},1,'last'), find(panelLoc{2},1,'first')               ...
                              ; cellfun(@sum,panelLoc)                                                ...
                              ];
        position            = [ gridX(panelSpan(1,2))         + axesMargin(3)                         ... left
                              , gridY(panelSpan(1,1))         + axesMargin(2)                         ... bottom
                              , panelStride(2)*panelSpan(2,2) - panelSpacing(2)-sum(axesMargin(3:4))  ... width
                              , panelStride(1)*panelSpan(2,1) - panelSpacing(1)-sum(axesMargin(1:2))  ... height
                              ];
        obj.gridIndex(iPanel) = sub2ind(size(panelGrid), panelSpan(1,1), panelSpan(1,2));

        %% Create panel
        obj.axs(iPanel)     = axes( 'Parent'                  , obj.fig                               ...
                                  , 'Units'                   , 'inches'                              ...
                                  , 'Position'                , position*PaneledFigure.PIX_TO_INCHES  ...
                                  , 'XColor'                  , FormatDefaults.axisCl                 ...
                                  , 'YColor'                  , FormatDefaults.axisCl                 ...
                                  , 'Layer'                   , 'top'                                 ...
                                  , 'ActivePositionProperty'  , 'Position'                            ...
                                  , 'FontSize'                , round(PaneledFigure.TICK_FONTSIZE * fontMultiplier)                       ...
                                  , 'LabelFontSizeMultiplier' , PaneledFigure.LABEL_FONTSIZE/PaneledFigure.TICK_FONTSIZE * fontMultiplier ...
                                  , 'TitleFontSizeMultiplier' , PaneledFigure.LABEL_FONTSIZE/PaneledFigure.TICK_FONTSIZE * fontMultiplier ...
                                  );
        if exist('enhanceCopying', 'file')
          enhanceCopying(obj.axs(iPanel));
        end
        
      end
      
      %% Panel grouping (controls whether the subsequent panels have a new label)
      obj.newGroup          = true(1,nPanels);
      obj.groupLabel        = gobjects(0);
      
    end
    
    %----- Draw hints for how to index panels; press a key/mouse button to continue
    function hint(obj)
      
      hHint             = gobjects(numel(obj.axs),2);
      for iPanel = 1:numel(obj.axs)
        if ~ishghandle(obj.axs(iPanel))
          continue;
        end
        
        xRange          = get(obj.axs(iPanel), 'XLim');
        yRange          = get(obj.axs(iPanel), 'YLim');
        hHint(iPanel,1) = text( obj.axs(iPanel), mean(xRange), mean(yRange)       ...
                              , sprintf('%d', iPanel)                             ...
                              , 'FontSize'    , 6*PaneledFigure.LABEL_FONTSIZE    ...
                              , 'FontWeight'  , 'bold'                            ...
                              , 'Color'       , [0 0 0]                           ...
                              , 'HorizontalAlignment'   , 'center'                ...
                              , 'VerticalAlignment'     , 'middle'                ...
                              );
        hHint(iPanel,2) = text( obj.axs(iPanel), mean(xRange), mean(yRange)       ...
                              , sprintf('%d', iPanel)                             ...
                              , 'FontSize'    , 5*PaneledFigure.LABEL_FONTSIZE    ...
                              , 'Color'       , [191 120 227]/255                 ...
                              , 'HorizontalAlignment'   , 'center'                ...
                              , 'VerticalAlignment'     , 'middle'                ...
                              );
      end
      
      waitforbuttonpress;
      delete(hHint);
      
    end
    
    %----- Export figure to disk, and individual panels to a sub-folder
    function export(obj, codeFile, info, hideEmptyPanels, mainPngOnly, labelOffset, varargin)
      
      %% Default arguments
      if nargin < 4 || isempty(hideEmptyPanels)
        hideEmptyPanels = true;
      end
      if nargin < 5 || isempty(mainPngOnly)
        mainPngOnly     = PaneledFigure.developmentMode();
      end
      if nargin < 6 || isempty(labelOffset)
        labelOffset     = zeros(size(obj.axs));
      end
      
      varargin          = [PaneledFigure.EXPORT_OPTIONS, varargin];
      labelOptions      = PaneledFigure.formatOptions('panelLabel');
      warnState         = warning('query', 'MATLAB:LargeImage');
      warning('off', 'MATLAB:LargeImage');
      
      %% Ensure exactly one label (vertical) offset value per panel
      if size(labelOffset,1) == 1
        labelOffset     = repmat(labelOffset, size(obj.panelGrid,1), 1);
      end
      if size(labelOffset,2) == 1
        labelOffset     = repmat(labelOffset, 1, size(obj.panelGrid,2));
      end
      
      %% Don't show panels with no contents
      if hideEmptyPanels
        for iPanel = 1:numel(obj.axs)
          if ishghandle(obj.axs(iPanel)) && isempty(get(obj.axs(iPanel), 'Children'))
            set(obj.axs(iPanel), 'Visible', 'off');
          end
        end
      end
      
      %% Add panel label per group
      delete(obj.groupLabel);
      obj.groupLabel    = gobjects(0);
      groupID           = obj.firstLabel;
      for iPanel = 1:numel(obj.axs)
        if ~ishghandle(obj.axs(iPanel)) || ~obj.newGroup(iPanel) || (hideEmptyPanels < 2 && strcmp(get(obj.axs(iPanel), 'Visible'), 'off'))
          continue;
        end
        
        %% Make space for y-ticks, if present
        if isempty(get(obj.axs(iPanel), 'YTick')) || strcmp(get(obj.axs(iPanel), 'YColor'), 'none')
          xOffset       = 0.335;
        else
          xOffset       = 0.65;
        end
        yText           = get(get(obj.axs(iPanel), 'YLabel'), 'String');
        if iscell(yText)
          xOffset       = xOffset + 0.3*(numel(yText) - 1);
        end
        
        %% Get axis dimensions in inches
        axisUnits       = get(obj.axs(iPanel), 'Units');
        set(obj.axs(iPanel), 'Units', 'Inches');
        axisPos         = get(obj.axs(iPanel), 'Position');
        set(obj.axs(iPanel), 'Units', axisUnits);
        
        %% Account for reversed axes in label placement
        iGrid           = obj.gridIndex(iPanel);
        xRange          = get(obj.axs(iPanel), 'XLim');
        yRange          = get(obj.axs(iPanel), 'YLim');
        if strcmp(get(obj.axs(iPanel), 'XDir'), 'reverse')
          labelX        = PaneledFigure.offsetLocation( xOffset, xRange, 2, axisPos(3), get(obj.axs(iPanel),'XScale'));
          xAlign        = 'right';
        else
          labelX        = PaneledFigure.offsetLocation(-xOffset, xRange, 1, axisPos(3), get(obj.axs(iPanel),'XScale'));
          xAlign        = 'left';
        end
        if strcmp(get(obj.axs(iPanel), 'YDir'), 'reverse')
          labelY        = PaneledFigure.offsetLocation(-labelOffset(iGrid), yRange, 1, axisPos(4), get(obj.axs(iPanel),'YScale'));
          yAlign        = 'top';
        else
          labelY        = PaneledFigure.offsetLocation( labelOffset(iGrid), yRange, 2, axisPos(4), get(obj.axs(iPanel),'YScale'));
          yAlign        = 'baseline';
        end
        
        %% Create and store label
        obj.groupLabel(end+1)                                                             ...
                        = text( double(labelX), double(labelY), char(groupID)             ...
                              , 'Parent'              , obj.axs(iPanel)                   ...
                              , 'FontSize'            , PaneledFigure.PANEL_FONTSIZE      ...
                              , 'HorizontalAlignment' , xAlign                            ...
                              , 'VerticalAlignment'   , yAlign                            ...
                              , 'Clipping'            , 'off'                             ...
                              , labelOptions{:}                                           ...
                              );
        groupID         = groupID + 1;
      end
      
      %% Store figure 
      outPath           = regexprep(codeFile, '[.][mM]$', '', 'once');
      outPrefix         = regexp(outPath, '^.*?[/\\]*([^/\\]*)?$', 'tokens', 'once');
      outPrefix         = outPrefix{:};
      if ~exist(outPath, 'dir')
        mkdir(outPath);
      end
      
      fprintf('  -->>  %s ...           ', outPath);
      set(obj.fig, 'UserData', info);
      savefig(obj.fig, fullfile(outPath, [outPrefix '.fig']), 'compact');
      
      
      %% Print info struct to a text file, if provided
      if ~isempty(info)
        fileID          = fopen(fullfile(outPath, [outPrefix '.txt']), 'w');
        PaneledFigure.printInfoStruct(fileID, info);
        fclose(fileID);
      end
      
      %% Store plot and individual groups of panels
      hasChanged        = nan(1, numel(obj.groupLabel)+1);
      for iOut = 1:numel(PaneledFigure.EXPORT_FORMAT)
        hasChanged(end) = PaneledFigure.exportIfDifferent ( obj.fig, [outPath, PaneledFigure.EXPORT_FORMAT{iOut}]   ...
                                                          , hasChanged(end), '-nocrop', varargin{:}                 ...
                                                          );

        %% Short-circuit in case of abbreviated output
        if mainPngOnly
          continue;
        end
      
        %% Loop over panels but store only per group
        panels          = gobjects(0);
        iGroup          = 0;
        for iPanel = 1:numel(obj.axs) + 1
          if iGroup <= numel(obj.groupLabel)
            fprintf('\b\b\b\b\b\b\b\b\b\b%4s %2d/%-2d', PaneledFigure.EXPORT_FORMAT{iOut}, iGroup, numel(obj.groupLabel));
            drawnow;
          end
          
          if iPanel > numel(obj.axs) || obj.newGroup(iPanel)
            %% Store the previous panel group
            if ~isempty(panels)
              label     = get(obj.groupLabel(iGroup), 'String');
              hasChanged(iGroup)                                                                                    ...
                        = PaneledFigure.exportIfDifferent ( panels, fullfile(outPath, ['panel_', label(1), PaneledFigure.EXPORT_FORMAT{iOut}])  ...
                                                          , hasChanged(iGroup), varargin{:}                         ...
                                                          );
              panels    = gobjects(0);
            end
            iGroup      = iGroup + 1;
          end
          
          %% Add panels to the current group
          if iGroup > numel(obj.groupLabel)
            break;
          end
          if iPanel <= numel(obj.axs)
            panels(end+1) = obj.axs(iPanel);
          end
        end
      end
      
      %% Restore warnings
      fprintf(' panels\n');
      
      warning(warnState.state, 'MATLAB:LargeImage');
      
    end
    
  end
    
  %________________________________________________________________________
  methods

    %----- Get method for panel, so that it is made the current axes; if a second output is requested, increments the index before returning the panel
    function [axs, index] = panel(obj, index)
      if nargout > 1
        index   = index + 1;
        while ~ishghandle(obj.axs(index))
          if index >= numel(obj.axs)
            error('PaneledFigure:panel', 'No more panels available to select.');
          end
          index = index + 1;
        end
      end
      
      if index > numel(obj.axs)
        error('PaneledFigure:panel', 'Panel index %d out of bounds, should be in the range [1,%d].', index, numel(obj.axs));
      end
      axs       = obj.axs(index);
      if ~ishghandle(axs)
        error('PaneledFigure:panel', 'Panel %d does not exist (has been deleted).', index);
      end
      
      set(0      , 'CurrentFigure', obj.fig);
      set(obj.fig, 'CurrentAxes'  , axs);
    end
    
    %----- Returns the total number of panels
    function count = numPanels(obj)
      count   = numel(obj.axs);
    end

    %----- Reshapes panels to have the desired relative widths (default equal) but in a row
    %      occupying the same extent as before. Panels can be deleted by specifying relativeSize as
    %      an array with less number of elements than the total number of redistributed panels.
    function distributeHorizontally(obj, index, varargin)
      obj.axs = PaneledFigure.distributePanels(obj.axs, 1, 3, index, varargin{:});
    end
    
    %----- Reshapes panels to have the desired relative heights (default equal) but in a column 
    %      occupying the same extent as before. Panels can be deleted by specifying relativeSize as
    %      an array with less number of elements than the total number of redistributed panels.
    function distributeVertically(obj, index, varargin)
      obj.axs = PaneledFigure.distributePanels(obj.axs, 2, 4, index, varargin{:});
    end
    
  end
  
  %________________________________________________________________________
  methods (Static)

    %----- Function to print a structure to a text file in a comprehensible way
    function printInfoStruct(file, info)

      %%
      infoFields    = fieldnames(info)';
      nameLength    = max(cellfun(@numel, infoFields));
      
      %%
      for field = infoFields
        fprintf(file, '  %s\n', field{:});
        fprintf(file, '%s\n', repmat('=', 1, nameLength+4));
        PaneledFigure.printInfo(file, info.(field{:}), '  ');
        fprintf(file, '\n');
      end
      
    end
    
    %----- Helper function to recursively print a structure 
    function printInfo(file, info, prefix)

      %% Define formatting method based on the type of data
      switch class(info)
        case 'table'
          %% Deduce display widths
          nameWidth     = max(cellfun(@numel, info.Properties.VariableNames));
          if isempty(info.Properties.RowNames)
            info.Properties.RowNames  = arrayfun(@(x) sprintf('%d',x), 1:size(info,1), 'UniformOutput', false);
          end
          preWidth      = max(cellfun(@numel, info.Properties.RowNames)) + 2;
          
          %% Print each row as a separate info block
          for iRow = 1:size(info,1)
            pre         = sprintf('%s %*s  ', prefix, preWidth, ['[' info.Properties.RowNames{iRow}, ']']);
            for iCol = 1:size(info,2)
              data      = info{iRow,iCol};
              if isempty(data)
                data    = '';
              elseif iscell(data) && numel(data) == 1
                data    = data{:};
              end
              
              %%
              fprintf(file, '%s %*s : ', pre, nameWidth, info.Properties.VariableNames{iCol});
              if ischar(data)
                data    = regexprep(data, '\n', ['$0' prefix repmat(' ',1,nameWidth+6)], 'lineanchors');
                data    = strtrim(data);
                fprintf(file, '%s\n', data);
              else
                data    = arrayfun(@(x) sprintf('%.*g', PaneledFigure.NUMBER_PRECISION, x), data, 'UniformOutput', false);
                fprintf(file, '%s\n', strjoin(data, '   '));
              end
              pre       = sprintf('%s %*s  ', prefix, preWidth, '');
            end
            fprintf(file, '%s\n', pre);
          end
          
        case 'struct'
          %%
          nameWidth     = max(cellfun(@numel, fieldnames(info))) + numel(prefix) + 4*(numel(info) > 1);
          for field = fieldnames(info)'
            for iObj = 1:numel(info)
              if numel(info) > 1
                postfix = sprintf('(%d)',iObj);
              else
                postfix = '';
              end
              
              name      = sprintf('%-*s', nameWidth, [prefix field{:} postfix]);
              if isstruct(info.(field{:}))
                PaneledFigure.printInfo(file, info.(field{:}), [prefix field{:} postfix '.']);
              else
                PaneledFigure.printInfo(file, info.(field{:}), [name ' : ']);
              end
            end
          end

        case 'char'
          fprintf(file, '%s%s\n', prefix, info);
          
        case 'logical'
          if numel(info) > PaneledFigure.MAX_DISP_ITEMS/2
            logicLabel  = {'F', 'T'};
          else
            logicLabel  = {'false', 'true'};
          end
          fprintf(file, '%s%s\n', prefix, strjoin(logicLabel(1+info), ' '));

        otherwise
          if isnumeric(info) && ~isvector(info) && max(size(info)) <= PaneledFigure.MAX_DISP_LENGTH
            %% Print tensors one slice (3rd and onwards dimension) at a time
            [info, infoSize, colWidth, tableWidth]   ...
                        = PaneledFigure.formatCells(num2cell(info));
            for iPage = 1:infoSize(3)
              if iPage > 1
                fprintf(file, '%s---[ %2d ]%s\n', prefix, iPage, repmat('-',1,tableWidth-8));
              end

              for iRow = 1:infoSize(1)
                fprintf(file, prefix);
                for iCol = 1:infoSize(2)
                  fprintf(file, '  %-*s ', colWidth(iCol), info{iRow,iCol,iPage});
                end
                prefix  = repmat(' ', 1, numel(prefix));
                fprintf(file, '\n');
              end
            end
          
          elseif isnumeric(info) && numel(info) <= PaneledFigure.MAX_DISP_ITEMS
            %% Other arrays
            rep         = arrayfun(@(x) sprintf('%.*g', PaneledFigure.NUMBER_PRECISION, x), info, 'UniformOutput', false);
            fprintf(file, '%s%s\n', prefix, strjoin(rep, '   '));
            
          elseif iscellstr(info) && numel(info) <= PaneledFigure.MAX_DISP_ITEMS
            fprintf(file, '%s%s\n', prefix, strjoin(info, ', '));
            
          elseif iscell(info) && all(size(info) <= PaneledFigure.MAX_DISP_ITEMS) && all(cellfun(@(x) isnumeric(x) || ischar(x), info(:)))
            %% Print table one page (3rd and onwards dimension) at a time
            [info, infoSize, colWidth, tableWidth]   ...
                        = PaneledFigure.formatCells(info);
            colStart    = {' ', '|'};
            for iPage = 1:infoSize(3)
              if iPage > 1
                fprintf(file, '%s---[ %2d ]%s\n', prefix, iPage, repmat('-',1,tableWidth-8));
              end
              
              for iRow = 1:infoSize(1)
                fprintf(file, prefix);
                for iCol = 1:infoSize(2)
                  fprintf(file, '%s %-*s ', colStart{1 + (iCol>1)}, colWidth(iCol), info{iRow,iCol,iPage});
                end
                prefix  = repmat(' ', 1, numel(prefix));
                fprintf(file, '\n');
              end
            end
            
          else
            dataSize    = sprintf('%dx', size(info));
            fprintf(file, '%s%s %s\n', prefix, dataSize(1:end-1), class(info));
            
          end
      end
      
    end
    
    %----- Export figure with change diff and lazy updates
    function hasChanged = exportIfDifferent(handle, targetFile, hasChanged, varargin)
      
      %% Check for changes if state is not specified
      tempFile        = [];
      if isnan(hasChanged)
        if exist(targetFile, 'file')
          %% Target already exists, make a temporary output and check for changes
          [~,~,ext]   = parsePath(targetFile);
          tempFile    = [tempname() ext];
          export_fig(handle, tempFile, varargin{:});
          testImage   = imread(tempFile);
          origImage   = imread(targetFile);
          hasChanged  = ~isequal(testImage, origImage);

        else
          %% Target doesn't exist, must have changed
          hasChanged  = true;
        end
      end
      
      %% Don't do anything if there are no changes
      if ~hasChanged
        if ~isempty(tempFile)
          delete(tempFile);
        end
        return;
      end
      
      %% Create target if not already done
      if isempty(tempFile)
        export_fig(handle, targetFile, varargin{:});
      else
        movefile(tempFile, targetFile);
      end
    end
    
    %----- Convenience function to locate the largest monitor
    function index = largestMonitor()
      monitors    = get(0,'monitor');
      screenArea  = prod( monitors(:,3:end), 2 );
      [~,index]   = max(screenArea);
    end
    
  end
  
end
