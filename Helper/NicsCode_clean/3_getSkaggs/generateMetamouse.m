function [metamouse] = generateMetamouse(argin1DFFthresholding, ...
    input2log, input3gaussianSmoothing, input4trialsAverageMap, ...
    input5numberShuffles, input6imagingFilename, input7thresholdSig, ...
    input8dimensions, input9binEdges, input10randomDimMethod, ...
    input11randomBinEdges, input12taskType, input13whichROIs, ...
    input14gaussianSmoothingToggle, input15LR_split)
% 
% Description: Calculates the mutual information for given imaging files
% using the getSkaggs analysis. it also allows the analysis to split 
% between left- and right-choice trials, showing results for each set of
% trials. 
%
% Sample Call:
% metamouse = generateMetamouse([11 4], 'noLog', 0, 'keepTrials', 10, {'E22_20170215_30p.modeling_NEW.mat', 'E39_20171012_40per.modeling.mat', 'E43_20170720_70per_userSetSD5minDur0.modeling.mat', 'E47_20170927_70per_userSetSD5minDur0.modeling.mat', 'E48_20170808_60per_userSetSD5minDur0.modeling.mat'}, 2, {'Evidence', 'Position'}, {[-20:20], [0:10:300]}, {'ViewAngle', 'Evidence'}, {[], []}, 'towers', 'all', 'both', 1);
% 
% *** INPUTS ***
% input1DFFthresholding - array of two integers, first integer sets the
% windowsize and second integer sets the thresholding
% if [0 0] == nothing happens
% if [0 4] == no filtering, threshold data 4*robustSTD
% if [11 4]== filter with windowsize=11 and threshold data 4*robustSTD
% if [11 0]== filter with windowsize=11, but don't threshold
%
% input2log - character vector
% 'log10' - take log10 of activity values with log10(v + 1);
% 'noLog' - do not take log10 activity
%
% input3gaussianSmoothing - double
% if ==0, no smoothing
% if >0, use value as sigma for gaussian smoothing
%
% input4trialsAverageMap - vector of integers, or character vector 'all'
% specify which trials should be used to construct the average maps for
% each neuron; useful for comparing correct vs error trials or left vs
% right trials. if specified as 'keepTrials', then construct average maps
% from only mainTrial and goodQuality trials from score.trial
%
%   ex:
%   'keepTrials' - do trials that are true for isGood and mainMaze
%   'keepTrials + leftChoice' - do 'keepTrials' and left choice trials
%   'keepTrials + rightChoice' - do 'keepTrials' and right choice trials
%
% input5numberShuffles - integer 0 or >0
% specify number of shuffles to retrieve/write, using circshift within
% trial
% if == 0, do not do shuffle tests
% if >0, do that number of shuffles
%
% input6imagingFilename - character vector filename for animal session
% example: 'E22_20170215_30p.modeling_NEW.mat'
%
% input7thresholdSig - threshold for significance of comparing skaggs
% values with shuffled skaggs values, e.g. 2 for 2SD above mean - ALSO
% used for the "shroud", i.e. threshold map, in pixelwise
%
% input8dimensions - cell array of character vectors specifying which
% dimension(s) to test for skaggs metric. case sensitive and must be found
% in output table of extractVariables()
% example: {'Evidence', 'Position', 'Choice'}
%
% input9binEdges - cell array of bin edges passed to discretize().
% ex: {[-20:20], [0:10:300]} for an ExY analysis
%
% input10randomDimMethod - 1-by-2 cell array of character vectors to
% specify which dimension to randomize with respect to a second dimension.
% Dim{1} is randomized with respect to its joint distribution with Dim{2}.
% ex: {'Evidence', 'Position'}  ...which is interpreted as 'randomly sample
% Evidence with respect to Position'
%
% input11randomBinEdges - 1-by-2 cell array of double vector. This sets the
% bin edges for the random dimensions in the same fashion of how
% input11binEdges sets the binEdges for the real dimension(s)
%
% input12taskType - 'towers' or 'alternation', used for extractVariables function
%
% input13whichROIs - which ROIs to use for Skaggs metric and average map
% calculations. Specify either 'all' or a double vector of ROI labels
%   'all' - use all ROIs in modeling.mat file
%   double vector - use only ROIs in vector
%
% input14gaussianSmoothingToggle - character vector to specify which
% part(s) of the analysis to use gaussian filtering on:
%       - 'skaggsOnly' allows smoothing on skaggs calculation and sequence
%       plot
%       - 'pixelwiseOnly' allows smoothing on bin-wise analysis of mean DFF
%       in pixelwise
%       - 'both' allows smoothing on both analyses above
%
% input15LR_split - split the analysis between L and R choice trials
%   0 - don't split
%   1 - do split; overrides input4trialsAverageMap and uses 'keepTrials +
%   leftChoice' and 'keepTrials + rightChoice'


imagingFiles = input6imagingFilename;

metamouse(length(imagingFiles)).animalID = [];
metamouse(length(imagingFiles)).filename = [];
if input15LR_split == 0
    metamouse(length(imagingFiles)).getSkaggs = [];
end
if input15LR_split == 1
    metamouse(length(imagingFiles)).LC = [];
    metamouse(length(imagingFiles)).RC = [];
end

strDims = [];
for i = 1:length(input8dimensions)
    if i > 1, strDims = [strDims, 'x']; end %#ok<*AGROW>
    switch input8dimensions{i}
        case 'Evidence', strDims = [strDims, 'E'];
        case 'Position', strDims = [strDims, 'Y'];
        case 'ViewAngle', strDims = [strDims, 'VA'];
        case 'Choice', strDims = [strDims, 'C'];
        otherwise, strDims = [strDims, input8dimensions{i}];
    end
end

if input15LR_split == 0
    for i = 1:length(imagingFiles)
        % Edit 2021/2/3 to be able to use fnameStruct
        if isstruct(imagingFiles)==1
            filename = imagingFiles(i).fname;
        else
            filename = imagingFiles{i};
        end
        animalID = strsplit(filename, '_');
        animalID = animalID{1};
        animalID = strsplit(animalID,'\');
        animalID = animalID{end};
        
        saveString = [animalID, '_', strDims,'_[', ...
            num2str(argin1DFFthresholding(1)), '_', num2str(argin1DFFthresholding(2)), ...
            ']_%d.mat'];
        
        if strcmp(input4trialsAverageMap, 'keepTrials + leftChoice')
            saveString = [animalID, '_', strDims,'_LC_[', ...
                num2str(argin1DFFthresholding(1)), '_', num2str(argin1DFFthresholding(2)), ']_%d.mat'];
        end
        if strcmp(input4trialsAverageMap, 'keepTrials + rightChoice')
            saveString = [animalID, '_', strDims,'_RC_[', ...
                num2str(argin1DFFthresholding(1)), '_', num2str(argin1DFFthresholding(2)), ']_%d.mat'];
        end
        
        metamouse(i).filename = filename;
        metamouse(i).animalID = animalID;
        
        instAnimal = ...
            getSkaggs(argin1DFFthresholding, input2log, input3gaussianSmoothing, ...
            input4trialsAverageMap, input5numberShuffles, saveString, ...
            1, filename, input7thresholdSig, input8dimensions, input9binEdges, ...
            input10randomDimMethod, input11randomBinEdges, input12taskType, ...
            input13whichROIs, input14gaussianSmoothingToggle);
        
        metamouse(i).getSkaggs = instAnimal;
    end
end
if input15LR_split == 1
    metamouseLC = generateMetamouse(argin1DFFthresholding, input2log, input3gaussianSmoothing, ...
        'keepTrials + leftChoice', input5numberShuffles, input6imagingFilename, input7thresholdSig,...
        input8dimensions, input9binEdges, input10randomDimMethod, input11randomBinEdges, ...
        input12taskType, input13whichROIs, input14gaussianSmoothingToggle, 0);
    
    metamouseRC = generateMetamouse(argin1DFFthresholding, input2log, input3gaussianSmoothing, ...
        'keepTrials + rightChoice', input5numberShuffles, input6imagingFilename, input7thresholdSig,...
        input8dimensions, input9binEdges, input10randomDimMethod, input11randomBinEdges, ...
        input12taskType, input13whichROIs, input14gaussianSmoothingToggle, 0);
    
    for i = 1:length(imagingFiles)
        metamouse(i).animalID = metamouseLC(i).animalID;
        metamouse(i).filename = metamouseLC(i).filename;
        metamouse(i).LC = metamouseLC(i).getSkaggs;
        metamouse(i).RC = metamouseRC(i).getSkaggs;
    end
end
