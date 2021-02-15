function outputBarSEM = nieh_barSEM(varargin)

counter=1;
for i=1:length(varargin)
    
    if size(varargin{i},1)>1 && size(varargin{i},2)>1
        for j=1:size(varargin{i},2)
            curData = varargin{i};
            curBar = curData(:,j);
            meanBar(counter) = nanmean(curBar);
            semBar(counter) = nieh_sem(curBar);
            counter = counter + 1;
        end
    else
        
        curBar = varargin{i};
        meanBar(counter) = nanmean(curBar);
        semBar(counter) = nieh_sem(curBar);
        counter = counter + 1;
        
    end
end


bar(meanBar, 'EdgeColor', 'none');
hold on;

er = errorbar([1:length(semBar)], meanBar, zeros(1,length(semBar)), semBar);

er.LineStyle = 'none';
er.CapSize = 0;

outputBarSEM.meanBar = meanBar;
outputBarSEM.semBar  = semBar;
outputBarSEM.ax1 = gca;

