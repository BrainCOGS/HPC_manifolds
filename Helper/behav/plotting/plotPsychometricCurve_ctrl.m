
function handle = plotPsychometricCurve_ctrl(perfPsych,plotWhat,ca,cl,plotFit,sh,applyDefaults)

if nargin < 2
  plotWhat = 'bin';
end
if nargin < 3 || isempty(ca)
  f  = figure;
  ca = gca(f);
end
if nargin < 4 || isempty(cl)
  cl = analysisParams.ctrlCl;
end
if nargin < 5
  plotFit = 1;
end
if nargin < 6 || isempty(sh)
  sh = analysisParams.ctrlShade;
end
if nargin < 7 || isempty(applyDefaults)
  applyDefaults = true;
end

try
  axis(ca);
catch
  figure(ca);
end

switch plotWhat
  case 'all'
    x    = perfPsych.perfPsych_xaxis;
    y    = perfPsych.perfPsychJ;
    l    = -perfPsych.perfPsychJSEM(:,1)+y;
    u    = perfPsych.perfPsychJSEM(:,2)-y;
    xl   = [-perfPsych.perfPsych_xaxis(end)-1 perfPsych.perfPsych_xaxis(end)+1];
    if plotFit
      fitx = perfPsych.fitAll.xaxis;
      fity = perfPsych.fitAll.curve;
    end
  case 'bin'
    x    = perfPsych.perfPsych_xaxisBins;
    y    = perfPsych.perfPsychJ_binned;
    l    = -perfPsych.perfPsychJ_binnedSEM(:,1)+y;
    u    = perfPsych.perfPsychJ_binnedSEM(:,2)-y;
    xl   = [-perfPsych.perfPsych_xaxisBins(end)-1 perfPsych.perfPsych_xaxisBins(end)+1];
    if plotFit(1) && ~isempty(perfPsych.fit)
      fitx = perfPsych.fit.xaxis;
      fity = perfPsych.fit.curve;
    else
      fitx = [];
      fity = [];
    end
end
y(x==0)=nan;
hold on
plot([0 0],[0 100],'--','color',[.7 .7 .7])
plot([x(1)-.5 x(end)+.5],[50 50],'--','color',[.7 .7 .7])
if plotFit(1)
%  handle = errorbar(x,y*100,l*100,u*100,'o','color',cl,'markersize',4,'markerfacecolor',cl);
  idx0 = x==0;
  x(idx0)=[]; y(idx0)=[]; l(idx0)=[]; u(idx0)=[];
  handle = errbar(x,y'*100,[l'*100;u'*100],cl,.75,0,'none');
  plot(x,y*100,'.','color',cl,'markersize',7);
  handle(end+1) = plot(fitx,fity*100,'-','color',sh,'linewidth',1);
else
  y(x==0) = [];
  l(x==0) = [];
  u(x==0) = [];
  x(x==0) = [];
  if numel(plotFit) > 1
    if plotFit(2)
      handle  = errorbar(x,y*100,l*100,u*100,'o-','color',cl,'markersize',4,'markerfacecolor',cl);
    else
      handle  = errorbar(x,y*100,l*100,u*100,'o','color',cl,'markersize',4,'markerfacecolor',cl);
    end
  else
    handle  = errorbar(x,y*100,l*100,u*100,'o-','color',cl,'markersize',4,'markerfacecolor',cl);
  end

end

box off;
xlim([x(1)-1 x(end)+1]);
ylim([0 100])

if applyDefaults
  set(gca,'fontsize',12,'xcolor','k','ycolor','k','ytick',0:25:100)
  ylabel('Went right (%)','fontsize',14,'color','k')
  xlabel('#R - #L towers','fontsize',14,'color','k')
end

