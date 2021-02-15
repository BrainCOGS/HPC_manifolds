
% This is a script that takes the .model.mat files of a animals and
%   1) fits N manifolds while leaving ONE random trial out
%   2) if needed, generates folder-structure to store the data
%
%  At its core, it uses the 'submitJobs' script that submits a matlab job
%  into the slurm engine.

taskType = 'towers';

dembed       = '2,3,4,5,6,7';   % embedding dimension
input_rng    = 1;               % random seed
Ncrossval    = 10;              % number of leave-out trials for crossvalidation
leaflist     = 500;             % list of leafs for fitting manifolds of different complexity
lm_fraclist  = 1;               % list of landmark_fractions

path              = '/jukebox/tank/enieh/mind/FINAL/Towers/';           % Storage path for the dF/F files
storage_directory = '/jukebox/tank/enieh/mind/FINAL/Towers/CrossValidation/';            % Path where data is to be stored


for animal = {'E48'}
    animal = char(animal);
    
    % Get the data
    if strcmp(animal, 'E39')
        raw_datafile = 'E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat';
    elseif strcmp(animal, 'E22')
        raw_datafile = 'E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat';
    elseif strcmp(animal, 'E43')
        raw_datafile = 'E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat';
    elseif strcmp(animal, 'E44')
        raw_datafile = 'E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat';
    elseif strcmp(animal, 'E47')
        raw_datafile = 'E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat';
    elseif strcmp(animal, 'E48')
        raw_datafile = 'E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat';
    elseif strcmp(animal, 'E65')
        raw_datafile = 'E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat';
    else
        error("no .modeling.mat file found")
    end

    % Making sure the folder and animal-subfolder for this class of problem exist
    if ~isfolder(storage_directory)
          disp(storage_directory);
          mkdir(storage_directory);
    end
    savepath = [storage_directory, animal '/'];
    if ~isfolder(savepath)
        disp(savepath);
        mkdir(savepath);
    end

    % General stuff for each animal: set rng, find trials, make folder if
    % needed
    fname = [path raw_datafile];
    rng(input_rng);
    nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);
    list_trials = unique(nic_output.behavioralVariables.Trial);
    rnd_idx = randperm(length(list_trials));
    random_order = list_trials(rnd_idx);
    
    folder_c = [savepath, 'allTrials/'];
    if ~isfolder(folder_c)
        disp(folder_c);
        mkdir(folder_c);
    end
    
    % And send it to queue for different random held-out trials
    for crossval_idx = 1:min([Ncrossval, length(list_trials)])
        leaveout_trial = random_order(crossval_idx);
        use_trials = random_order;
        use_trials(crossval_idx) = [];
        folder_t = [folder_c, 'trialout_', num2str(leaveout_trial), '/'];
        if ~isfolder(folder_t)
            disp(folder_t);
            mkdir(folder_t);
        end
        
        %%% and send it into the queue for different nleafs
        files = dir(folder_t);
        for nleafs = leaflist
            for lmf = lm_fraclist
                %test whether this has already been calculated
                searchstr = [num2str(nleafs), '_lmf_', num2str(lmf), '_manifold.mat'];
                job_already_done = length(strfind([files.name], searchstr))>0;
                if job_already_done
                    disp('Job already done. Nothing to do here.')
                else
                    disp([folder_t, ' --- ', num2str(nleafs)]);
                    disp('... was missing - job submitted!')
                    use_trials_sorted = sort(use_trials');
                    use_trials_str = regexprep(num2str(use_trials_sorted),'\s+',',');
                    trial_predict = num2str(leaveout_trial);
                    input_minLeaves = num2str(nleafs);
                    command_string = "mindExecutable('" + ...
                                use_trials_str      + ...
                        "','" + trial_predict       + ...
                        "','" + dembed              + ...
                        "','" + num2str(input_rng)  + ...
                        "','" + input_minLeaves     + ...
                        "','" + num2str(lmf)        + ...
                        "','" + fname               + ...
                        "','" + folder_t            + ...
                        "','" + taskType + "')";

                    submitJobs(command_string, taskType);
                end
            end
        end
    end
end


    
