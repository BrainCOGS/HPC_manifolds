%%
%
function [centers] = toBinCenters(edges)

  if numel(edges) < 2
    error('At least two bin edges must be provided.');
  end
  
  centers = ( edges(1:end-1) + edges(2:end) ) / 2;
  
end

