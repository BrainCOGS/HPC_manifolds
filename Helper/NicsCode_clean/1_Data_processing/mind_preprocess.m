function clean_dataDFF = mind_preprocess(dataDFF, toggleSmoothing, toggleThreshold)
% windowSize of 11 gives a 733 ms window
% if ~toggleSmoothing==0, do the smoothing using the value as the window
% size
% if ~toggleThreshold==0, do the thresholding, using the value as the
% number of robustSTD to set the threshold to

% Add gaussian smoothing to the data preprocessing
if ~toggleSmoothing==0
    dataDFF = smoothdata(dataDFF, 'gaussian',toggleSmoothing);
end

% Threshold dataDFF values
if ~toggleThreshold==0
    for i=1:size(dataDFF,2)
        dataDFFROI = dataDFF(:,i);
        threshold = robustSTD(dataDFFROI)*toggleThreshold;
        dataDFFROI(dataDFFROI<threshold)=0;
        dataDFF(:,i) = dataDFFROI;
    end
end

clean_dataDFF = dataDFF;
