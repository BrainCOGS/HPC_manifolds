function chosenROIs = mind_findTopROIs(fname, chosenI, topI, rngInput)

rng(rngInput)

nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);
ROIactivities = nic_output.ROIactivities;
sumActivity = sum(ROIactivities);
[~, i2] = sort(sumActivity, 'descend');
topIs = i2(1:topI);
chosenROIs = topIs(randperm(topI,chosenI));

figure;
for i=1:chosenI
    subplot(chosenI,1,i)
    plot(ROIactivities(:,chosenROIs(i)));
end
