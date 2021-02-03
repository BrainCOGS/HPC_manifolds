function handle = errbar(x,y,err,cl,linew,horizFlag,linestyle)

% handle = errbar(x,y,err,cl,linew,horizFlag,linestyle)

if nargin < 3; error('I need at least 3 inputs');       end
if nargin < 4 || isempty(cl); cl = 'k';                 end
if nargin < 5 || isempty(linew); linew = .75;           end
if nargin < 6 || isempty(horizFlag); horizFlag = false; end
if nargin < 7 || isempty(linestyle); linestyle = '-';   end

if size(x,1) > size(x,2)
  x = x';
end
if size(y,1) > size(y,2)
  y = y';
end
if size(err,1) > size(err,2)
  err = err';
end
if size(err,1) == 1
  err = [err; err];
end

hold on

if horizFlag
  plot(x,y,'o','color',cl,'markerfacecolor',cl); 
  handle = plot([x-err(1,:); x+err(2,:)],[y; y],'-','linewidth',linew,'color',cl); 
else
  plot(x,y,'linestyle',linestyle,'linewidth',linew,'color',cl); 
  handle = plot([x; x],[y-err(1,:); y+err(2,:)],'-','linewidth',linew,'color',cl); 
end