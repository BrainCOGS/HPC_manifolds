%--------------------------------------------------------------------------
%%
function [y] = fiteval(f, x, nonzeroOnly)

  %-----  Input parsing
  if nargin < 3
    nonzeroOnly = false;
  end
  

  %-----  Assume all structures are C++ splines
  if isstruct(f)
    % Special case with an array of f corresponding to each x
    if numel(f) > 1 && numel(f) == numel(x)
      y         = zeros(size(x));
      for iX = 1:numel(x)
        y(iX)   = evaluateSpline(x(iX), f(iX), nonzeroOnly);
      end
    else
      y         = evaluateSpline(x, f, nonzeroOnly);
    end
    
  %-----  Otherwise assume it is Matlab handle or cfit object
  else
    if nonzeroOnly
      select    = ( x ~= 0 );
    else
      select    = true(size(x));
    end

    y           = zeros(size(x));
    if iscell(f) && numel(f) > 1 && numel(f) == numel(x)
      for iX = 1:numel(x)
        if select(iX)
          y(iX) = arrayfun(f(iX), x(iX));
        end
      end
    else
      y(select) = f(x(select));
    end
  end
  
end
