function dat = fitManifold2Video(fname, videopath, savepath, savepath_manifold, do_collectrawdata)
    % Fits a manifold to the raw pixel values of the video
    %
    % Requires as data
    % - XXX.modeling.mat of the session
    % - Videos of all the individual trials (pull from dj)
    % - Details of behavior. This can be skipped if done once (also pulled
    % form dj)
    % 
    %
    % Usage:
    %    fname:             Path to .modeling.mat of that session
    %    videopath:         Path to folder with all the videofiles
    %    savepath:          Path where to find/store the raw_data (frames & behavior but no manifold). Computation requires datajoint
    %    savepath_manifold: Path where the manifold will be saved
    %    do_collectrawdata: 
    %                       If yes: 
    %                           raw_data is collected and saved at `savepath`
    %                       If no:
    %                           `savepath` is loaded
    %
    % Saves:
    %  - the manifold from the video in the usual struct.
    %
    % Example use:
    %         fitManifold2Video('/Users/ms81/project_GeometryCognition/HPC_data/E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat', ...
    %             '/Users/ms81/E65_trials/', ...
    %             '/Users/ms81/TT_embedding_rawdata.mat', ...
    %             '/Users/ms81/TT_embedding_manifold.mat', ...
    %             true)
    %

    rng(42);
        
    if do_collectrawdata
        
        % Get data from dj
        setenv('DJ_HOST', 'datajoint00.pni.princeton.edu')
        dj.conn()
        setenv('DB_PREFIX', 'u19_')
        key.subject_fullname = 'hnieh_E65';
        key.session_date = '2018-02-02';
        relevant_data = behavior.TowersBlockTrial & key;
        [blocks, n_trials] = fetchn(behavior.TowersBlock & key, 'block', 'n_trials');
        cumtrialsblock = cumsum([0;n_trials]);

        % Get the data from the modeling.mat file:
        out = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 11, [11 4], fname,'none','towers', 1, 1);
        use_trials = unique(out.behavioralVariables.Trial);    
        load(fname)  %to get score
        
        sampling = 17;  %very 17th pixel to make computation faster
        training_set = use_trials(1:round(length(use_trials)*0.95));
        
        all_data = [];
        all_test_data = [];
        data_positionY = [];
        data_positionX = [];
        data_leftT     = [];
        data_rightT    = [];
        data_luminance = [];
        data_phase     = [];
        data_viewangle = [];
        data_trial     = [];
        data_evidence  = [];

        for tt = use_trials'
            
            block = score.trial(tt).block;
            trial = score.trial(tt).index - cumtrialsblock(block);
            
            file = ['hnieh_E65-2018-02-02-B', + num2str(block), '-T', num2str(trial), '.mp4'];
            disp(['Reading... ' file])

            %open the file, and get all the frames
            v = VideoReader([videopath, file]);
            maxframes = floor( v.Duration * v.FrameRate )-1;
            data_video = zeros(maxframes, 63*105);  %63x105 from scaled fullHD
            luminance = zeros(maxframes, 1);

            for idx = 1:maxframes
                frame = read(v,idx);
                luminance(idx) = mean(mean(frame(:,:,3)));
                [sx, sy, ~] = size(frame);
                dx = int16(1:sampling:sx-sampling);
                dy = int16(1:sampling:sy-sampling);
                av_frame = zeros(length(dx), length(dy));
                for av = 0:sampling-1
                    av_frame = av_frame + double(frame(dx+av,dy+av,3));
                end
                B = av_frame/(sampling-1);  % 3 because we are interested in the "B" channel
                data_video(idx,:) = B(:);
            end

            if sum(tt==training_set)>0
            % attached to the collection
                all_data = [all_data; data_video];
            else
                all_test_data = [all_test_data; data_video];
            end
            
            %Get the data from the database
            data_thistrial = relevant_data & ['block =', num2str(block)] & ['trial_idx =', num2str(trial)];

            %Reorganize so that it can be plotted
            [onL, offL, onR, offR, position, iterations] = fetchn(data_thistrial, 'cue_onset_left', 'cue_offset_left', 'cue_onset_right', 'cue_offset_right', 'position', 'iterations' );
            onL = onL{1};
            offL = offL{1};
            onR = onR{1};
            offR = offR{1};
            evi = zeros(int16(v.Duration * v.FrameRate), 1);

            L = zeros(int16(v.Duration * v.FrameRate), 1);
            for idx =  1:length(onL)
                if offL(idx) < iterations  %make sure this actually happened
                    L(onL(idx) : offL(idx)) = 1;
                    evi(onL(idx):end) = evi(onL(idx):end) + 1;
                end
            end
            R = zeros(int16(v.Duration * v.FrameRate), 1);
            for idx =  1:length(onR)
                if offR(idx) < iterations
                    R(onR(idx) : offR(idx)) = 1;
                    evi(onR(idx):end) = evi(onR(idx):end) - 1;
                end
            end

            pos = position{1};

            data_viewangle  = [data_viewangle;  pos(1:maxframes,3)];
            data_positionY  = [data_positionY;  pos(1:maxframes,2)];
            data_positionX  = [data_positionX;  pos(1:maxframes,1)];

            data_leftT     = [data_leftT;     L(1:maxframes)];
            data_rightT    = [data_rightT;    R(1:maxframes)];
            data_luminance = [data_luminance; luminance(1:maxframes)];
            data_phase     = [data_phase;     (1:maxframes)'];
            data_trial     = [data_trial;     ones(size((1:maxframes)'))*tt];
            data_evidence  = [data_evidence;  evi ];
            
            [i_cue_entry, i_mem_entry, i_turn_entry, iterations] = fetchn(data_thistrial, 'i_cue_entry', 'i_mem_entry', 'i_turn_entry', 'iterations');
            Tblock = score.trial(tt).block;
            Ttrial = tt - min(find([score.trial.block] == block))+1;
            
            %Make health check
            if score.trial(tt).iterations ~= iterations  % Test 1-4: is dj and score consistent?
                disp("break")
            elseif score.trial(tt).iTurnEntry ~= i_turn_entry
                disp("break")
            elseif score.trial(tt).iMemEntry  ~= i_mem_entry
                disp("break")                
            elseif score.trial(tt).iCueEntry  ~= i_cue_entry
                disp("break")     
            elseif abs(iterations - maxframes) > 5          % Test 5-6: is dj and video consistent?
                disp("break")
            elseif abs(length(luminance) - length(pos)) > 5
                disp("break")
            elseif Tblock ~= block                          % Test 7-8 are trial/block consistent?
                disp("break")
            elseif Ttrial ~= trial
                disp("break")
            end
            
        end
        save(savepath,'-v7.3')
    end
   

        
    %% Fit Manifold
    load(savepath)
    traind = ismember(data_trial, training_set);
    flag = (data_phase(traind)<450) & (data_positionY(traind)>0);
    data_for_mind = all_data(flag, :); %Something like "17k Timepoints x 7k neurons"
    [timepoints, ~] = size(data_for_mind);
    
    mindparameters.dt = 1;
    mindparameters.pca.n = 0.95;
    mindparameters.dim_criterion = .95;
    mindparameters.ndir = 2;
    mindparameters.min_leaf_pts = 500;
    mindparameters.ntrees = 100;
    mindparameters.verbose = true;
    mindparameters.lm.n = 1000;

    mindparameters.rwd.type = 'discrete';
    mindparameters.rwd.sym = 'avg';
    mindparameters.rwd.all_geo = true;
    mindparameters.rwd.d = 2;
    mindparameters.rwd.var_scale = 0.1;

    mindparameters.embed.type = 'rwe';
    mindparameters.embed.d = nan;
    mindparameters.embed.mode = 'mds';
    mindparameters.embed.local = false;

    mindparameters.learnmapping = true;
    mindparameters.mapping.k = [1:10, 15:5:50];
    mindparameters.mapping.lambda = [0, 10.^(-8:.5:0)];
    mindparameters.mapping.mode = 'lle';
    mindparameters.mapping.nfolds_lle = 10;

    mindparameters.prune_lm_by_time = false;

    dembed = [1,2,3,4,5,6,7]; % embedding dimensions

    data = struct();
    data.t = reshape(1:timepoints,timepoints,1);
    data.f = data_for_mind;

    result = struct();
    result.forestdat = mindAsFunction(data, mindparameters);
    result.mindparameters = mindparameters;

    embedparameters = mindparameters; 
    embedparameters.embed.d = dembed;
    [~, result.allembed] = embedAsFunction(result.forestdat, embedparameters);

    [t,n] = size(result.forestdat.pca.f);
    disp("Needed " + num2str(n) + " PCs for 95% variance")

    save(savepath_manifold,'-v7.3')
    
    %%% Plot manifold
    % pca_coords = result.forestdat.pca.model.transform(data_for_mind, 0.95);
    % y = result.allembed(3).f2m.map.transform(pca_coords);
    % scatter3(y(:,1), y(:,2), y(:,3), [], data_luminance(flag), '.')
    
end