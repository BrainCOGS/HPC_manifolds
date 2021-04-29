function outputLuminance2HPC = luminance2HPC(fnameStruct, videopath, vfilename)

% Get behavioral data for E65 from Edwards files

load(fnameStruct(7).fname) %to get score
nic_output = extractVariables('all', 2, 'keepTrials', 2, 0, 0, 5, [11 4], fnameStruct(7).fname, 'none', 'towers', 1, 1);
trialn = unique(nic_output.behavioralVariables.Trial)'; 

%% 
trial_luminance = [];
trial_position  = [];
for tt = trialn
        
    data_thistrial_ed = nic_output.behavioralVariables(nic_output.behavioralVariables.Trial == tt,:);
    
    block = score.trial(tt).block;
    trial = tt - min(find([score.trial.block] == block))+1;
    pos  = score.trial(tt).position;
      
    %Open the file, and measure luminance
    file = [vfilename, num2str(block), '-T', num2str(trial), '.mp4'];
    v = VideoReader([videopath, file]);
    maxframes = floor( v.Duration * v.FrameRate )-1;
    luminance = zeros(maxframes, 1);
    for idx = 1:maxframes    %in steps of five, a bit less data
        frame = read(v,idx);
        luminance(idx) = mean(mean(frame(:,:,3))); %mean luminance of blue channel
    end
    
    for idx = 1:length(data_thistrial_ed.Position)
        y_here = data_thistrial_ed.Position(idx);

        [~, minidx] = min((pos(:,2) - y_here).^2);
        trial_luminance = [trial_luminance, luminance(minidx)];
        
    end
   disp(tt)
end

outputLuminance2HPC.luminance = luminance;
outputLuminance2HPC.trial_luminance = trial_luminance;

