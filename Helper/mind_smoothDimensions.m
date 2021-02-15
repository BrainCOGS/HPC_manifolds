function newData = mind_smoothDimensions(oldData, behaviorData, numSmooth)

for i=1:size(oldData,1)
    
    distances = sqrt(sum(bsxfun(@minus, oldData, oldData(i,:)).^2,2));
    [~, indmin] = sort(distances,'ascend');
    closeData = behaviorData(indmin(1:numSmooth));
    newData(i) = sum(closeData)/length(closeData);
    
end