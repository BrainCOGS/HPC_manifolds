function submitJobs(to_execute_input, taskType, animalID)
    % sends job to the cluster
    % to_execute should be a STRING which is the full matlab command
    
    to_execute = to_execute_input{1}(1:end-2) + "','$SLURM_JOB_ID')";  % Append slurmid to job
    
    [~, user]= system("echo `whoami`");
    
    % resources on spock:
    % https://npcdocs.princeton.edu/index.php/SpockResources
    % Max 24 cores, 250GB, max. 12GB/core
    
    if contains(user, "ms81")
        JobPath = '/usr/people/ms81/logs/';
%         This works for most jobs, but not E22/39:
         partition = 'Brody';
         MaxRuntime = '2880';
         MaxMemory = '40000';
         CpuTask   =  '3';
    else
        if exist('animalID')
            if strcmp(animalID,'E22')==1 || strcmp(animalID,'E39')==1 || strcmp(animalID,'E48')==1
                JobPath = '/jukebox/tank/enieh/mind/logs/';
                partition = 'Brody';
                MaxRuntime = '5760';
                MaxMemory = '120000';
                CpuTask   =  '8';
            else
                JobPath = '/jukebox/tank/enieh/mind/logs/';
                partition = 'Brody';
                MaxRuntime = '5760';
                MaxMemory = '30000';
                CpuTask   =  '2';
            end
        else
            JobPath = '/jukebox/tank/enieh/mind/logs/';
            partition = 'all';
            MaxRuntime = '5750';
            MaxMemory = '120000';
            CpuTask   =  '8';
        end
    end
    
    strsp = split(to_execute,"',");
    heldout = strsp{2};
    rng = strsp{3};
    minLeaves = strsp{5};
    lmf = strsp{6};
    
    if strcmp(taskType,'towers')==1 || strcmp(taskType,'Towers')==1
        strsp = split(to_execute,"Towers/");
    elseif strcmp(taskType,'T7')==1
        strsp = split(to_execute,"T7/");
    elseif strcmp(taskType,'Alternation')==1
        strsp = split(to_execute,"Alternation/");
    elseif strcmp(taskType,'AlternationJeff')==1
        strsp = split(to_execute,"AlternationJeff/");
    end
    
    animalname = strsp{2};
    
    if length(heldout) < 10
        infos = [animalname(1:3), '_HO_', heldout(2:end), '_ML_', minLeaves(2:end), '_LF_', lmf(2:end)];
    else
        infos = [animalname(1:3), '_HO_', 'ALL', '_ML_', minLeaves(2:end), '_LF_', lmf(2:end)];
    end

    JobName = ['MIND_%j_', infos];
    LogFiletype = '.out';

    Header = sprintf([ 'date\n',                                     ...
            'echo "In the directory: `pwd` "\n',                     ...
            'echo "As the user: `whoami` "\n',                       ...
            'echo "on host: `hostname` "\n',                         ...
            'echo "With access to cpu id(s): "\n',                        ...
            'cat /proc/$$/status | grep Cpus_allowed_list\n']);
    ProgramModule = 'module load matlab/R2018a';
    Program = 'matlab -nosplash -nojvm -nodisplay -nodesktop -r ';
    Shell = '#!/usr/bin/env bash';
    coderepo = "addpath(genpath('./mind'))";
    Argument = sprintf('"try; %s; %s; ', coderepo, to_execute) + "catch me; fprintf('%s / %s\n',me.identifier,me.message); end; exit";

    %Add time and day to jobname and logfile name

    LogfileName = sprintf('%s%s%s',JobPath,JobName,LogFiletype);
    disp(LogfileName);
    %add shell
    RunCommand = sprintf( "%s\n\n", Shell);
    %add parameters
    RunCommand = sprintf( "%s#SBATCH -J %s \n#SBATCH -o %s \n#SBATCH -p %s \n#SBATCH -t %s \n#SBATCH --mem %s \n#SBATCH --cpus-per-task %s",...
        RunCommand, JobName, LogfileName, partition, MaxRuntime, MaxMemory, CpuTask);
    %add program
    RunCommand = sprintf("%s\n\n%s", RunCommand, Header);
    RunCommand = sprintf( '%s\n\n%s\n%s%s"', RunCommand, ProgramModule, Program, Argument);
    
    %This executes the command via slurm
    fid = fopen('tmp.sh','wt');
    fprintf(fid, '%s', RunCommand);
    fclose(fid);
    system("sbatch ./tmp.sh");
    pause(1)
    system("rm ./tmp.sh");
    
end




                              



