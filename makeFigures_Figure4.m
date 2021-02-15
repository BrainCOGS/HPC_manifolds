
%% Generate the doublets and triplets

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Either run this code
% Or load:

outputDoublets_E22 = findDoublets(11,fnameStruct(1).fname,[5 0], 1);
outputDoublets_E39 = findDoublets(11,fnameStruct(2).fname,[5 0], 1);
outputDoublets_E43 = findDoublets(5, fnameStruct(3).fname,[5 0], 1);
outputDoublets_E44 = findDoublets(5, fnameStruct(4).fname,[5 0], 1);
outputDoublets_E47 = findDoublets(5, fnameStruct(5).fname,[5 0], 1);
outputDoublets_E48 = findDoublets(5, fnameStruct(6).fname,[5 0], 1);
outputDoublets_E65 = findDoublets(5, fnameStruct(7).fname,[5 0], 1);

outputTriplets_E22 = findTriplets(outputDoublets_E22);
outputTriplets_E39 = findTriplets(outputDoublets_E39);
outputTriplets_E43 = findTriplets(outputDoublets_E43);
outputTriplets_E44 = findTriplets(outputDoublets_E44);
outputTriplets_E47 = findTriplets(outputDoublets_E47);
outputTriplets_E48 = findTriplets(outputDoublets_E48);
outputTriplets_E65 = findTriplets(outputDoublets_E65);

%% Plot the example doublets (Fig. 4a)

load(fnameStruct(2).fname);
deltaT = score.deltaT;
plotSequenceTrialROI(219,outputDoublets_E39, fnameStruct(2).fname, deltaT);
plotSequenceTrialROI(122,outputDoublets_E39, fnameStruct(2).fname, deltaT);


%% Make the metamouse and generate plots (Fig. 4b-c, f-j 

makeOutputDoublets_Metamouse
makeOutputTriplets_Metamouse

plotDoublets(outputDoublets_metamouse);
plotTriplets(outputTriplets_metamouse);

rightSig_triplet(:,2) = [outputTriplets_metamouse.saveAll_triplets_Sig_Right.cell1and2only_predict];
rightSig_triplet(:,1)   = [outputTriplets_metamouse.saveAll_triplets_Sig_Right.prediction];
rightSig_triplet(:,3) = [outputTriplets_metamouse.saveAll_triplets_Sig_Right.cell3only_predict];
rightSig_triplet(any(isnan(rightSig_triplet), 2), :) = [];
[~,pr1] = ttest(rightSig_triplet(:,1), rightSig_triplet(:,2))
[~,pr2] = ttest(rightSig_triplet(:,1), rightSig_triplet(:,3))
[~,pr3] = ttest(rightSig_triplet(:,2), rightSig_triplet(:,3))

leftSig_triplet(:,2) = [outputTriplets_metamouse.saveAll_triplets_Sig_Left.cell1and2only_predict];
leftSig_triplet(:,1)   = [outputTriplets_metamouse.saveAll_triplets_Sig_Left.prediction];
leftSig_triplet(:,3) = [outputTriplets_metamouse.saveAll_triplets_Sig_Left.cell3only_predict];
leftSig_triplet(any(isnan(leftSig_triplet), 2), :) = [];
[~,pl1] = ttest(leftSig_triplet(:,1), leftSig_triplet(:,2))
[~,pl2] = ttest(leftSig_triplet(:,1), leftSig_triplet(:,3))
[~,pl3] = ttest(leftSig_triplet(:,2), leftSig_triplet(:,3))

[p, pr] = ttest([outputTriplet_metamouse.saveAll_triplets_Sig_Left.cell1and2only_predict], [outputTriplet_metamouse.saveAll_triplets_Sig_Left.trip1_2_no3_shuffle])
[p, pr] = ttest([outputTriplet_metamouse.saveAll_triplets_Sig_Right.cell1and2only_predict], [outputTriplet_metamouse.saveAll_triplets_Sig_Right.trip1_2_no3_shuffle])


%% To plot the trajectory between two cells in doublets (Fig. 4d)

outputTrajectory_E39 = mind_plotDoubletTrajectory(fnameStruct(2).fname_mani, fnameStruct(2).fname, 1, 11, [5 0], 3, outputDoublets_E39, 219, [6 7], 0);

%% Make scatter plot of distance between doublet events (Fig. 4e)

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

% Either regenerate doublets or load

% Either run this code
% Or load: 

outputDoubletManiCorrs_E22 = calcDoubletManifoldCorrelation(fnameStruct(1).fname, fnameStruct(1).fname_mani, 100, 5, outputDoublets_E22);
outputDoubletManiCorrs_E39 = calcDoubletManifoldCorrelation(fnameStruct(2).fname, fnameStruct(2).fname_mani, 100, 5, outputDoublets_E39);
outputDoubletManiCorrs_E43 = calcDoubletManifoldCorrelation(fnameStruct(3).fname, fnameStruct(3).fname_mani, 100, 5, outputDoublets_E43);
outputDoubletManiCorrs_E44 = calcDoubletManifoldCorrelation(fnameStruct(4).fname, fnameStruct(4).fname_mani, 100, 5, outputDoublets_E44);
outputDoubletManiCorrs_E47 = calcDoubletManifoldCorrelation(fnameStruct(5).fname, fnameStruct(5).fname_mani, 100, 5, outputDoublets_E47);
outputDoubletManiCorrs_E48 = calcDoubletManifoldCorrelation(fnameStruct(6).fname, fnameStruct(6).fname_mani, 100, 5, outputDoublets_E48);
outputDoubletManiCorrs_E65 = calcDoubletManifoldCorrelation(fnameStruct(7).fname, fnameStruct(7).fname_mani, 100, 5, outputDoublets_E65);

