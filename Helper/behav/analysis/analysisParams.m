classdef analysisParams
    
    properties (Constant)
        
        % paths
        savepath                = '/Users/lucas/Documents/Princeton/data/laserGalvo/';
        savepathBehav           = '/Users/lucas/Documents/Princeton/data/behav/';
        gridpath                = '/Users/lucas/Documents/Princeton/code/laserCtrl/grid/';
        serverRoot              = '/Volumes/braininit/RigData/';
        serverPathBehav         = '/Volumes/braininit/RigData/VRLaser/behav/lucas/';
        serverPathLsr           = '/Volumes/braininit/RigData/VRLaser/LaserGalvo1/';
        pathForSpock            = '/jukebox/braininit/RigData/';
        pathForSpockBehav       = '/jukebox/braininit/RigData/VRLaser/behav/lucas/';
        pathForSpockLsr         = '/jukebox/braininit/RigData/VRLaser/LaserGalvo1/';
        serverPathLsrAnalysis   = '/Volumes/braininit/Analysis/laserGalvo/';
        serverPathLsrAnalysisPC = '\\bucket.pni.princeton.edu\braininit\Analysis\laserGalvo\';
        pathForSpockLsrAnalysis = '/jukebox/braininit/Analysis/laserGalvo/';
        serverPathBehavModels   = '/Volumes/braininit/Analysis/BehavModels/';
        cpFileTypes             = {'notes','cnsl','im'}; % file types to copy from server unprocessed
        
        % list of all mice
        mice            = { 'wt12';'wt13';'wt11';'wt6';'vg1';'vg2';'vg3';'vg4';'vg5';       ...
                            'vg6';'vg7';'vg8';'vg9';'vg10';'vg11';'vg12';'vg14';'vg15';     ...
                            'vg16';'vg17';'vg18';'vg20';'vg22';'vg23';'vg24';'vg25';        ...
                            'vg26';'vg27';'vg28';'vg29';'vg30';'vg31';'vg32';'vg33';'vg34'; ...
                            'vg36';'vg37';'vg38';'vg41';'vg42';'vg39';'vg40';'vg51';'vg53'
                          }; %;'wt5';
        
        % code used in logs
        rightCode       = 1
        leftCode        = 0
        nilCode         = -1
        abortCode       = NaN
        dotsPerRev      = 2605.9 % laser rig dots calibration
        
        % analysis parameters
        perfTh          = .6; % default performance cutoff
        psychbins       = -15:5:15; % -15:2:15;%psychometric curve x axis
        bootalpha       = .025; % alpha level for bootstrap (pre-correction)
        
        % plotting 
        lsrCl               = [57 133 255]./255;
        ctrlCl              = [0 0 0];
        RlsrCl              = [0 0 1];
        RctrlCl             = [0 0 .5];
        LlsrCl              = [1 0 0];
        LctrlCl             = [.5 0 0];
        lsrShade            = [151 193 252]./255;
        ctrlShade           = [.7 .7 .7];
        matchTrialCl        = [73 133 27]./255;
        matchTrialShade     = [86 184 144]./255;
        matchLsrTrialCl     = [138 43 226]./255;
        matchLsrTrialShade  = [153 80 224]./255;
        mazeColors          = {[.75 .75 .75],[.5 .5 .5],[.25 .25 .25],[0 0 0],[1 .7 .2],[.6 0 0],[.8 0 0],[.8 .3 .3],[1 0 0],[.9 .2 .9],[.5 .2 .9],[.2 .5 1],[.2 .2 1],[0 0 1]};
        mazeColorsCond      = {[.75 .75 .75],[.5 .5 .5],[.25 .25 .25],[0 0 0],[.6 0 0],[1 .7 .2],[1 0 0],[.9 .2 .9],[.5 .2 .9],[0 0 .4],[.6 .6 1]};
        colormap            = 'red2blue'%'hot'
        cbarLim             = [-1e-5 1e-5];
        cbarLimEs           = [-.5 .5];
        radiusMM            = [.08 1.6]; 
        ctSkullPath         = '/Users/lucas/Documents/Princeton/figs&movies/cartoonTransparentSkull.tiff'; 
        ctBrainPath         = 'allenBrainOverlay.tiff';%'/Users/lucas/Documents/Princeton/figs&movies/allenBrainOverlay.tiff'; 
        ctSkullPxlPerMM     = 30.5; 
        ctBrainPxlPerMM     = 41.4; 
        ctSkullBregma       = [146 143]; 
        ctBrainBregma       = [203 274];
        myOrange            = [245 145 32]./255;
        myBlue              = [57 133 255]./255;
        myRed               = [255 0 0]./255;
        cueHalfCl           = [.3 .3 .3];%[188 19  254]./255;
        cueHalfShade        = [220 150 255]./255;
        testcl              = [0 0 .1];
        
        % for behavioral analysis
        behavProtType   = {'PoissonBlocksReboot3m';'PoissonBlocksCondensed3m';'PoissonBlocksCondensed3m_Ben';'PoissonBlocksReboot3mTransientB';'PoissonBlocksReboot3mTransient'};
        miceBehav       = {'vg1';'vg2';'vg3';'vg4';'vg5';'vg6';'vg7';'vg8';'vg9';'vg10';'vg11';'vg12';'vg14';'vg15';'vg17';...
                            'wt6';'wt11';'wt12';'wt13';'ty2';...
                            'ai1';'ai2';'ai3';'ai4';'ai5';'ai6';...
                            'k39';'k40';'k41';'k42';'k43';'k44';'k45';'k46';...
                            'k49';'k50';'k51';'k52';'k53';'k54';'k55';'k56';'k57';'k58';...
                            'B490';'B441';'B706';'B737';'B705';'B440';...
                            'ai7';'ai8'};
        miceBehavLegacy = {'vg1';'vg2';'vg3';'vg4';'vg5';'vg6';'vg7';'vg8';'vg9';'vg10';'vg11';'vg12';'vg14';'vg15';...
                            'wt6';'wt11';'wt12';'wt13';...
                            'ai1';'ai2';'ai3';'ai4';'ai5';'ai6';...
                            'k39';'k40';'k41';'k42';'k43';'k44';'k45';'k46';...
                            'k49';'k50';'k51';'k52';'k53';'k54';'k55';'k56';'k57';'k58';...
                            };
        genotypes       = {'vgat';'wt';'ai93';'thy1';'datcre';'yfp'};
        genCl           = {[0 0 1],[.3 .3 .3],[0 .5 0],[.7 0 .5]};
        areaCl          = [[0 129 59]./255;[0 224 205]./255;[0 116 225]./255;    ...
                           [231 118 22]./255;[231 179 42]./255;[172 51 59]./255; ...
                           [207  79 223]./255;[63 140 109]./255;[221 114 144]./255]
        multiCl         = { ...
                            [  .750      .750      .750], ...
                            [  .500      .500      .500], ...
                            [  .250      .250      .250], ...
                            [     0         0         0], ...
                            [0.2540         0         0], ...
                            [0.5079         0         0], ...
                            [0.7619         0         0], ...
                            [1.0000    0.0159         0], ...
                            [1.0000    0.2698         0], ...
                            [1.0000    0.5238         0], ...
                            [1.0000    0.7778         0], ...
                            [1.0000    1.0000    0.0476], ...
                            [1.0000    1.0000    0.4286], ...
                            [1.0000    1.0000    0.8095], ...
                            [1.0000    1.0000    0.9365]};                 
                        
    end
    
    %%
    properties
        
        filters
        
    end
    
    %%
    methods
      
      function applyAxisDefaults(obj,axisHandle,cl)
        if nargin < 3; cl = 'k'; end
        set(axisHandle,'ycolor',cl,'xcolor',cl,'zcolor',cl,'fontsize',12,'box','off')
      end
      
      function applyAxisLbls(obj,axisHandle,xlbl,ylbl,titlestr,zlbl)
        if nargin < 5; titlestr = []; end
        if nargin < 6; zlbl     = []; end
        axes(axisHandle)
        xlabel(xlbl,'fontsize',14)
        ylabel(ylbl,'fontsize',14)
        if ~isempty(zlbl)     
          zlabel(zlbl,'fontsize',14);                                    
        end
        if ~isempty(titlestr) 
          title(titlestr,'fontsize',15,'fontweight','bold','color',get(axisHandle,'xcolor')); 
        end
      end
      
      function applyFigDefaults(obj,figHandle,panels,cl)
        if nargin < 4; cl = 'w'; end
        set(figHandle,'color',cl,'position',[10 10 min([1200 800; panels.*250])])
      end
      
      function rootdir = getRootDir(obj,spockFlag,localFlag)
        if nargin < 3; localFlag = false; end
        if spockFlag
          rootdir   = analysisParams.pathForSpockLsrAnalysis;
        else
          if ispc
            rootdir = analysisParams.serverPathLsrAnalysisPC;
          else
            if localFlag
              rootdir = analysisParams.savepath;
            else
              rootdir = analysisParams.serverPathLsrAnalysis;
            end
          end
        end
      end
      
    end
end
