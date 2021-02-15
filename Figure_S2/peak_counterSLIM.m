function [num_peaks, pixelwiseInput] = peak_counterSLIM(pixelwiseInput, minInput)

min_pixels_for_peak = minInput;

num_peaks = zeros(1,length(pixelwiseInput));
for i = 1:length(pixelwiseInput)
    data = pixelwiseInput(i).realSignifMap;
    CC = bwconncomp(data>0);
    for j = 1:CC.NumObjects
        if length(CC.PixelIdxList{j}) > min_pixels_for_peak
            num_peaks(i) = num_peaks(i) + 1;
        else
            data(CC.PixelIdxList{j})=0;
        end
    end
    
    pixelwiseInput(i).realSigniMap_minusSmall = data;
end

figure
histogram(num_peaks, [-0.5:1:6.5], 'Normalization','probability')
xlabel('number of peaks')
ylabel('frequency')
set(gca,'box','off')