function outputDecodePCAVariable = mind_decodePCAVariable(fnameStruct, taskType, varTypeList)

% varTypeList is two vars that will be tested

%% Track Inputs and Set up RNG and config variables

outputDecodePCAVariable.argins.fnameStruct = fnameStruct;
outputDecodePCAVariable.argins.taskType    = taskType;
outputDecodePCAVariable.argins.varTypeList = varTypeList;


rngInput = 42;
rng(rngInput);

dimEmbed = 5;
numFolds = 5;

outputRegressOutViewAngle2.config.rngInput = rngInput;
outputRegressOutViewAngle2.config.dimEmbed = dimEmbed;
outputRegressOutViewAngle2.config.numFolds = numFolds;


%%
all_coefs = [];
for i=1:length(fnameStruct)
    
    disp(i)
    load(fnameStruct(i).fname_mani);
    fname = fnameStruct(i).fname;
    
    if strcmp(taskType, 'alternation')==1 || strcmp(taskType, 'Alternation')==1
        nic_output = extractVariables('all', 2, 'goodTrials', 2, 0, 0, 11, [11 4], fname,'none','alternation', 1, 1);
        
    elseif strcmp(taskType, 'towers')==1 || strcmp(taskType, 'tower')==1 || strcmp(taskType, 'Towers')==1 || strcmp(taskType, 'T7')==1
        nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none', 'towers', 1, 1);
        
    elseif strcmp(taskType, 'alternationJeff')==1 || strcmp(taskType, 'AlternationJeff')==1
        nic_output = extractVariables('all', 6, 'goodTrials', 2, 0, 0, 11, [11 4], fname,'none','alternation', 1, 1);
        
    end
    
    behavioralVariables = nic_output.behavioralVariables;
    behavioralVariables = behavioralVariables(outMind.Datarange,:);
    
    fn = fieldnames(behavioralVariables);
    var1 = table2array(behavioralVariables(:,strcmp(fn,varTypeList(1))));
    var2 = table2array(behavioralVariables(:,strcmp(fn,varTypeList(2))));
    
    % Get manifold
    sample = outMind.dataDFF;
    pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
    manifold3d = outMind.dat.allembed([outMind.dat.allembed.d]==dimEmbed).f2m.map.transform(pca_coords);
    
    % Change into orthogonal coordinates.
    [coeff,score,latent] = pca([var1, var2]);
    
    %And do the decoding of both
    numFolds = 5;
    len = floor(length(var2)/numFolds);
    coefs = zeros(2, numFolds);
    
    CV = generateCrossValSet_v2(behavioralVariables, numFolds);
    
    for cv = 1:numFolds  %crossvalidation block
        disp(['Fold... ', num2str(cv)])
        
        training_set = CV(cv).trainLocations;
        test_set     = CV(cv).testLocations;

        for t = [1,2]
            variables_train = score(training_set, t);
            variables_test  = score(test_set, t);

            Mdl = fitrgp(manifold3d(training_set, :), variables_train, 'Standardize', false);
            yfit = predict(Mdl, manifold3d(test_set, :));
            cc = corrcoef(variables_test, yfit);
            coefs(t,cv) = cc(2,1);
        end
    end

    disp(mean(coefs,2)')
    ani_coefs{i} = coefs;
    all_coefs = [all_coefs; mean(coefs,2)'];  % Collect the data
end

outputDecodePCAVariable.ani_coefs = ani_coefs;
outputDecodePCAVariable.all_coefs = all_coefs;

disp("Reconstruction, and p-values:")
for score_idx = [1,2]
    disp(['r=', num2str(mean(all_coefs(:,score_idx))), ' +- ', ...
    num2str(std(all_coefs(:,score_idx))/sqrt(length(all_coefs(:,score_idx)))) ]);

    disp(['p=', num2str(signrank(all_coefs(:,score_idx))) ] )  % signrank bc. correlations are not gaussian
    disp('---')
end

