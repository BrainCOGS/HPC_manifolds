function outputCompareTowersAlternationVA = mind_compareTowersAlternationVA(fnameStruct, fnameStruct_Alternation)

%% Calculate decoding of VA and Y in Alternation Task

for i=1:length(fnameStruct_Alternation)
    outputNonlinearDecoding_VA_ALT = mind_nonlinearDecoding_dimX_All(fnameStruct_Alternation(i).fname, fnameStruct_Alternation(i).fname_mani,5,'GP','ViewAngle','Alternation',0,1,[], 5);
    outputNonlinearDecodingAll_VA_ALT{i} = outputNonlinearDecoding_VA_ALT;
    meancorrAll_VA_ALT(i) = outputNonlinearDecoding_VA_ALT.meancorr;
    disp(['Animal ' num2str(i) ' of ' num2str(length(fnameStruct_Alternation)) ' finished']);
end

%% Calculate decoding of VA and Y in Towers Task

for i=1:length(fnameStruct)
    outputNonlinearDecoding_VA = mind_nonlinearDecoding_dimX_All(fnameStruct(i).fname, fnameStruct(i).fname_mani,5,'GP','ViewAngle','Towers',0,1,[], 5);
    outputNonlinearDecodingAll_VA{i} = outputNonlinearDecoding_VA;
    meancorrAll_VA(i) = outputNonlinearDecoding_VA.meancorr;
    disp(['Animal ' num2str(i) ' of ' num2str(length(fnameStruct)) ' finished']);
end

ranksum(meancorrAll_VA, meancorrAll_VA_ALT)

outputCompareTowersAlternationVA.meancorrAll_VA     = meancorrAll_VA;
outputCompareTowersAlternationVA.meancorrAll_VA_ALT = meancorrAll_VA_ALT;

