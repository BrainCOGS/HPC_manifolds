
%% Collects results

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

%% Setup
animalList = {'E22', 'E39', 'E43', 'E44', 'E47', 'E48', 'E65'};

config.input_rng = 1;           % random seed
Ncrossval        = 10;          % number of leave-out trials for crossvalidation
nleafs           = 500;         % list of leafs for fitting manifolds of different complexity
lmf              = 1;           % list of landmark_fractions

storage_directory = 'M:\enieh\mind\FINAL\Towers\CrossValidation\';

%% Main Code

counter = 1;
for animal = animalList
    
    animal = char(animal);
    
    % Get the data
    if strcmp(animal, 'E22')
        fname = fnameStruct(1).fname;
    elseif strcmp(animal, 'E39')
        fname = fnameStruct(2).fname;
    elseif strcmp(animal, 'E43')
        fname = fnameStruct(3).fname;
    elseif strcmp(animal, 'E44')
        fname = fnameStruct(4).fname;
    elseif strcmp(animal, 'E47')
        fname = fnameStruct(5).fname;
    elseif strcmp(animal, 'E48')
        fname = fnameStruct(6).fname;
    elseif strcmp(animal, 'E65')
        fname = fnameStruct(7).fname;
    else
        error("no .modeling.mat file found")
    end
    
    savepath = [storage_directory, animal '/'];
    
    rng(config.input_rng);
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);
    list_trials = unique(nic_output.behavioralVariables.Trial);
    rnd_idx = randperm(length(list_trials));
    random_order = list_trials(rnd_idx);
    
    folder_c = [savepath, 'allTrials/'];
    
    M = NaN(Ncrossval, 6);
    for crossval_idx = 1:min([Ncrossval, length(list_trials)])
        
        leaveout_trial = random_order(crossval_idx);
        folder_t = [folder_c, 'trialout_', num2str(leaveout_trial), '/'];
        
        %load the three embeddings;
        f = [folder_t, 'minLeaves_' num2str(nleafs), '_lmf_', num2str(lmf), '_crossval.mat'];
        disp(f);
        if isfile(f)
            load(f);
            M(crossval_idx,:) = output_crossVal.corrTrial;
        end
    end
    
    M2 = squeeze(mean(M,1));
    maxReconstruct(:,counter) = M2;
    counter = counter + 1;

end

