function outputNonlinearDecoding_VarGroup = mind_nonlinearDecoding_VarGroup(fnameStruct, numFolds, regType, varTypeList, taskType, shuffleToggle, manifoldToggle, dimEmbed)

% Uses mind_nonlinearDecoding_dimX_All
% varTypeList is a list of variables to be decoded, i.e. {'Evidence', 'Position'}

%% Set up config and argins
config.rng        = 42;
rng(config.rng);

argins.fnameStruct = fnameStruct;
argins.numFolds   = numFolds;
argins.regType    = regType;
argins.varTypeList    = varTypeList;
argins.taskType   = taskType;
argins.shuffleToggle = shuffleToggle;
argins.manifoldToggle = manifoldToggle;
argins.dimEmbed = dimEmbed;

outputNonlinearDecoding_VarGroup.config    = config;
outputNonlinearDecoding_VarGroup.argins    = argins;

%% Do the decoding

for j=1:length(varTypeList)
    
    varType = varTypeList{j};
    
    for i=1:length(fnameStruct)
        
        outputNonlinearDecoding = mind_nonlinearDecoding_dimX_All(fnameStruct(i).fname, fnameStruct(i).fname_mani,numFolds,regType,varType,taskType,shuffleToggle,manifoldToggle,[], dimEmbed);
        
        corrList(i,:) = outputNonlinearDecoding.corr2;
        corrListmean(i) = outputNonlinearDecoding.meancorr;
        
        outputNonlinearDecodingSave{i} = outputNonlinearDecoding;
        
        if strcmp(varType, 'PriorChoice') || strcmp(varType, 'PriorCorrect') || strcmp(varType, 'Choice') || strcmp(varType, 'ChoiceCorrect')
            meancorPerList(i) = outputNonlinearDecoding.meancorPer;
            binValmean_CVList(i,:) = outputNonlinearDecoding.binValmean_CV;
        end
        
    end
    
    outputNonlinearDecoding_VarGroup(j).corrList = corrList;
    outputNonlinearDecoding_VarGroup(j).corrListmean = corrListmean;
    
    if strcmp(varType, 'PriorChoice') || strcmp(varType, 'PriorCorrect') || strcmp(varType, 'Choice') || strcmp(varType, 'ChoiceCorrect')
        outputNonlinearDecoding_VarGroup(j).meancorPerList = meancorPerList;
        outputNonlinearDecoding_VarGroup(j).binValmean_CVList = binValmean_CVList;
    end
    
    outputNonlinearDecoding_VarGroup(j).outputNonlinearDecodingSave = outputNonlinearDecodingSave;
end
