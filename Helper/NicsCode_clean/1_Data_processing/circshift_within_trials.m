function [out] = circshift_within_trials(x, trials, shuffleIndLimit)
% Description: Shuffle method for data x. This function does a circular 
% shift on the DFF values of each ROI once per trial independently of all
% other ROIs. 
% 
% Sample Call: 
% out = circshift_within_trials(x, trials, 10);
% 
% ***** INPUTS *****
% x - M-by-N matrix with M frames of observations of N ROIs. 
% trials - M-by-1 vector that lists the trial ID of each M frame. 
% shuffleIndLimit - double; sets the range of indices to choose for
%   circshift as [1+shuffleIndLimit, end-shuffleIndLimit].
%
% ***** OUTPUTS *****
% out - x, with the circshift-by-trial method

out = NaN(size(x)); % preallocate output
uniqT = unique(trials); % list of trial IDs

for i = 1:length(uniqT) % loop through trials
    for j = 1:size(x, 2) % loop through ROIs
        % for each ROI, get the DFF values where the frame's trial is the
        % same as the i-index
        v = x(trials == uniqT(i), j);
        
        % find the shift index using shuffleIndLimit
        shiftInd = randi([1+shuffleIndLimit, ...
            length(v)-shuffleIndLimit], 1);
        
        % circshift v with the shift index
        v = circshift(v, shiftInd);
        
        % save the circshifted values in the output
        out(trials == uniqT(i), j) = v;
    end
end
end