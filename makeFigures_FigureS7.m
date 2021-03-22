
%% Plot the multiple doublet examples with heat maps

% First load the doublet from E39
load('C:\Neuroscience\imaging\FINAL\sequences_Data\outputDoubletsTriplets_all.mat')

% Then load the position mutual information data
load('C:\Neuroscience\imaging\FINAL\getSkaggs_Data\out_Y_all.mat')

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

randDoublet = [6 100 102 184 237 258 313 383 504 742 996 1113 1388 1499 1698 1862 2286 2539 2967 3359 3971 5295 5638 6249 6625];
plotSequenceTrialROI_multiple(outputDoublets_E39, fnameStruct(2).fname, out_E39_Y, randDoublet);