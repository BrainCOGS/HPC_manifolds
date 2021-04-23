function outputBarSEM = nieh_barSEMpaired(varargin)

for i=1:length(varargin)
   
    curBar = varargin{i};
    allpoints(:,i) = curBar;
    meanBar(i) = nanmean(curBar);
    semBar(i) = nieh_sem(curBar);
       
end

bar(meanBar, 'EdgeColor', 'none');
hold on;
er = errorbar([1:length(varargin)], meanBar, zeros(1,length(varargin)), semBar);
er.LineStyle = 'none'; 
er.CapSize = 0;

outputBarSEM.meanBar = meanBar;
outputBarSEM.semBar  = semBar;
outputBarSEM.ax1 = gca;

for i=1:size(allpoints,1)
plot(allpoints(i,:),'.-k');
end

