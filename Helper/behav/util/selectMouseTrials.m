% Convenience function to apply a mouseID selection criterion to either a subset or the full log.
% Example usage:
%   newLog              = selectMouseTrials(lg, mouseID);
%   [trialType, choice] = selectMouseTrials(lg, mouseID, 'trialType', 'choice');
% You can also ask for the selection vector as an additional output after the above:
%   [..., sel]          = selectMouseTrials(lg, mouseID, ...);
%
% The first form returns a struct in the same format as lg but with only the subset of trials
% matching the given mouse ID. 
%
% In the second form, a list of fields of interest is explicitly specified and only those data are
% returned post selection.
%
function varargout = selectMouseTrials(lg, mouseID, varargin)
  
  %% Special case for metamouse (for speed)
  if numel(mouseID) == 1 && mouseID < 0
    if isempty(varargin)
      varargout           = {lg};
    else
      varargout           = cell(size(varargin));
      for iVar = 1:numel(varargin)
        varargout{iVar}   = lg.(varargin{iVar});
      end
    end
    varargout{end+1}      = true(size(lg.trialType));
    return;
  end
  
  
  %% Apply mouseID trial selection
  sel                     = ismember(lg.mouseID, mouseID);
  
  if isempty(varargin)
    %% Copy entire log
    subset                = struct();
    for var = fieldnames(lg)'
      if numel(lg.(var{:})) == numel(lg.trialType)
        subset.(var{:})   = lg.(var{:})(sel);
      else
        subset.(var{:})   = lg.(var{:});
      end
    end
    varargout             = {subset, sel};
    
  else
    %% Copy a selected set of variables
    varargout             = cell(size(varargin));
    for iVar = 1:numel(varargin)
      if numel(lg.(varargin{iVar})) == numel(lg.trialType)
        varargout{iVar}   = lg.(varargin{iVar})(sel);
      else
        varargout{iVar}   = lg.(varargin{iVar});
      end
    end
    varargout{end+1}      = sel;
  end
  
end
