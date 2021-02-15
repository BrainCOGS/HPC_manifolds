function outputHPC2HPC = map_HPC2HPC(fnameStruct, bestfit_position, bestfit_evidence)

%% Set up some variables
t0 = clock;
target_dim = 5;


%% Analysis
M_position = zeros(length(fnameStruct));
M_evidence = zeros(length(fnameStruct));
for animal_i = 1:length(fnameStruct)
    
    fname = fnameStruct(animal_i).fname;
    fname_mani = fnameStruct(animal_i).fname_mani;
    
    data_i = get_data(fname, fname_mani, target_dim);
    models = make_Mdl(data_i);         % Fit the model to animal_i

    for animal_j = 1:length(fnameStruct)
        
        if animal_i == animal_j
            M_position(animal_i, animal_j) = bestfit_position(animal_i);
            M_evidence(animal_i, animal_j) = bestfit_evidence(animal_i);
        else
            
            rng(42); %rng here, so that we can easily reproduce every loop
            
            fname = fnameStruct(animal_j).fname;
            fname_mani = fnameStruct(animal_j).fname_mani;
            
            data_j = get_data(fname, fname_mani, target_dim);
            [data_train, data_test] = split_data(data_j);
            
            rotation_parameters = rotate_manifold(models, target_dim, data_train);
            fit = evaluate_manifold(models, target_dim, data_test, rotation_parameters);
            M_position(animal_i, animal_j) = fit.position_corrcoef;
            M_evidence(animal_i, animal_j) = fit.evidence_corrcoef;
        end
        
        animal_j
    end
    animal_i
end

duration = clock - t0;
disp(duration);

outputHPC2HPC.M_position = M_position;
outputHPC2HPC.M_evidence = M_evidence;

end

%% Helpers

function [data_train, data_test] = split_data(data)

    CV = generateCrossValSet_v2(data.behavioralVariables, 2);  % two folds
    trainset = ismember(data.trial, CV(1).train);
    testset  = ismember(data.trial, CV(1).test);
    
    data_train.manifold = data.all_manifold(trainset, :);
    data_train.position = data.position(trainset, :);
    data_train.evidence = data.evidence(trainset, :);

    data_test.manifold = data.all_manifold(testset, :);
    data_test.position = data.position(testset, :);
    data_test.evidence = data.evidence(testset, :);
end
            
function fit = evaluate_manifold(models, target_dim, data, x)

    Mdlposition = models.position;
    Mdlevidence = models.evidence;
    test_position = data.position;
    test_evidence = data.evidence;
    test_manifold = data.manifold;
    if target_dim == 5
        T = dim5rotations(x);
    elseif target_dim == 3
        T = dim3rotations(x);
    end
    posfit = predict(Mdlposition, (T*test_manifold')');
    evifit = predict(Mdlevidence, (T*test_manifold')');
    coef = corrcoef(posfit(isfinite(test_position)), test_position(isfinite(test_position)));
    fit.position_corrcoef = coef(2,1);

    coef = corrcoef(evifit(isfinite(test_evidence)), test_evidence(isfinite(test_evidence)));
    fit.evidence_corrcoef = coef(2,1);
end

function x = rotate_manifold(models, target_dim, data_train)

    Mdlposition = models.position;
    Mdlevidence = models.evidence;
    train_manifold = data_train.manifold;
    train_position = data_train.position;
    train_evidence = data_train.evidence;
    
    scale_estimate  = prctile(models.evidence.X(:), 90) ./ prctile(train_manifold(:), 90);
    if target_dim == 5
        x0 = [zeros(1,10), scale_estimate];
    elseif target_dim == 3
        x0 = [zeros(1,3), scale_estimate];
    end
    bb = @(x) minfu(x, Mdlposition, Mdlevidence, target_dim, train_manifold, train_position, train_evidence);
    options = optimset('MaxIter', 200, 'Display', 'iter', 'PlotFcns', @optimplotfval);
    x = fminunc(bb, x0, options);
end

function data = get_data(fname, fname_mani, target_dim)

    load(fname_mani)
    sample = outMind.dataDFF;
    pca_coords = outMind.dat.forestdat.pca.model.transform(sample, outMind.config_input.mindparameters.pca.n);
    data.all_manifold = outMind.dat.allembed([outMind.dat.allembed.d]==target_dim).f2m.map.transform(pca_coords);
    
    out = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname, 'none', 'towers', 1, 1);
    data.behavioralVariables = out.behavioralVariables;
    data.position = out.behavioralVariables.Position(outMind.Datarange,:);  % We want position
    data.evidence = out.behavioralVariables.Evidence(outMind.Datarange,:);  % We want evidence
    data.trial    = out.behavioralVariables.Trial(outMind.Datarange,:);     % For Crossval, we need trials

end

function models = make_Mdl(data)
    rng(42) % make this reproducible
    Mdlposition = fitrgp(data.all_manifold, data.position, 'Standardize', false);
    Mdlevidence = fitrgp(data.all_manifold, data.evidence, 'Standardize', false);
    models.position = Mdlposition;
    models.evidence = Mdlevidence;
    
end

function gof = minfu(x, Mdlposition, Mdlevidence, target_dim, all_manifold, position, evidence)

    if target_dim == 5
        T = dim5rotations(x);
    elseif target_dim == 3
        T = dim3rotations(x);
    end
    posfit = predict(Mdlposition, (T*all_manifold')');
    evifit = predict(Mdlevidence, (T*all_manifold')');
    coef = corrcoef(posfit(isfinite(position)), position(isfinite(position)));
    a = coef(2,1);
    
    coef = corrcoef(evifit(isfinite(evidence)), evidence(isfinite(evidence)));
    b = coef(2,1); 
    
    gof = -a-b;
end

function T = dim5rotations(x)
    % Calculate the 5-dim rotations from the special orthogonal group
    % With SO(n), one can rotate axis 1 into axes 2, 3,...,N. 
    % One axis 2, you can rotate it into 3,...N. 
    % So, SO(n) has (N-1)+(N-2)+...+1 = N(N-1)/2 generators. 
    % -> SO(2) has 1, 
    % -> SO(3) has 3, 
    % -> SO(4) has 6,
    % -> SO(5) has 10. 
    %
    % So for SO(5), I'm calculating the exponential map from its 
    % Lie algebra.
    
    % tests: det(T) is always 1,
    
    M1 = [0, -1, 0, 0, 0; 1, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
    M2 = [0, 0, 1, 0, 0; 0, 0, 0, 0, 0; -1, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
    M3 = [0, 0, 0, 0, 0; 0, 0, -1, 0, 0; 0, 1, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
    M4 = [0, 0, 0, -1, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0; 0, 0, 0, 0, 0];
    M5 = [0, 0, 0, 0, 0; 0, 0, 0, -1, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 0, 0, 0, 0];
    M6 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, -1, 0; 0, 0, 1, 0, 0; 0, 0, 0, 0, 0];
    M7 = [0, 0, 0, 0, -1; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0];
    M8 = [0, 0, 0, 0, 0; 0, 0, 0, 0, -1; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0];
    M9 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, -1; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0];
    M10 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, -1; 0, 0, 0, 1, 0];
    
    scale = x(11)*diag([1, 1, 1, 1, 1]);
%     scale = diag([x(11), x(12), x(13), x(14), x(15)]);

    T = expm(x(1)*M1)*expm(x(2)*M2)*expm(x(3)*M3)*expm(x(4)*M4)*expm(x(5)*M5)* ...
        expm(x(6)*M6)*expm(x(7)*M7)*expm(x(8)*M8)*expm(x(9)*M9)*expm(x(10)*M10)*scale;
end

function T = dim3rotations(x)
    % Same idea as with with dim=5
    
    M1 = [0, -1, 0; 1, 0, 0; 0, 0, 0];
    M2 = [0, 0, 1; 0, 0, 0; -1, 0, 0];
    M3 = [0, 0, 0; 0, 0, -1; 0, 1, 0];
    
    scale = diag([x(4), x(5), x(6)] );
    
    T = expm(x(1)*M1)*expm(x(2)*M2)*expm(x(3)*M3)*scale;
end