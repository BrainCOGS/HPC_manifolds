function fnameStruct = mind_makeFnameStruct(user1, taskType, cpuType)

% cpuType can be "laptop" or "spock"

if strcmp(user1,'Edward')==1
    
    if strcmp(taskType, 'towers')==1
        
        if strcmp(cpuType,'laptop') ==1
            
            fnameStruct(1).fname      = 'C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat';
            % fnameStruct(1).fname_mani = 'C:\Neuroscience\imaging\E22_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(1).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E22_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(1).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170227.mat';
            
            fnameStruct(2).fname      = 'C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat';
            % fnameStruct(2).fname_mani = 'C:\Neuroscience\imaging\E39_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(2).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E39_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(2).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171103.mat';
            
            fnameStruct(3).fname      = 'C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat';
            % fnameStruct(3).fname_mani = 'C:\Neuroscience\imaging\E43_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(3).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E43_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(3).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170802.mat';
            
            fnameStruct(4).fname      = 'C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat';
            % fnameStruct(4).fname_mani = 'C:\Neuroscience\imaging\E44_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(4).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E44_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(4).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171018.mat';
            
            fnameStruct(5).fname      = 'C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat';
            % fnameStruct(5).fname_mani = 'C:\Neuroscience\imaging\E47_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(5).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E47_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(5).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20170927.mat';
            
            fnameStruct(6).fname      = 'C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat';
            % fnameStruct(6).fname_mani = 'C:\Neuroscience\imaging\E48_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(6).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E48_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(6).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170829.mat';
            
            fnameStruct(7).fname      = 'C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat';
            % fnameStruct(7).fname_mani = 'C:\Neuroscience\imaging\E65_minLeaves_500_lmf_1_dim_2to7_manifold.mat';
            fnameStruct(7).fname_mani = 'C:\Neuroscience\imaging\FINAL\Towers\E65_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            fnameStruct(7).fname_behav = 'C:\Neuroscience\imaging\FINAL\Towers\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180202.mat';
            
        elseif strcmp(cpuType,'spock')==1
                 
            fnameStruct(1).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat';
            fnameStruct(1).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E22_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            
            fnameStruct(2).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat';
            fnameStruct(2).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E39_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            
            fnameStruct(3).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat';
            fnameStruct(3).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E43_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            
            fnameStruct(4).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat';
            fnameStruct(4).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E44_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            
            fnameStruct(5).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat';
            fnameStruct(5).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E47_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            
            fnameStruct(6).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat';
            fnameStruct(6).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E48_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
            
            fnameStruct(7).fname      = '/jukebox/tank/enieh/mind/FINAL/Towers/E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat';
            fnameStruct(7).fname_mani = '/jukebox/tank/enieh/mind/FINAL/Towers/E65_minLeaves_500_lmf_1_manifold_Towers_FINAL.mat';
        
        end
    elseif strcmp(taskType, 'T7')==1
        fnameStruct(1).fname      = 'C:\Neuroscience\imaging\FINAL\E22_20170125_30per_userSetSD11minDur0.modelingFINAL.mat';
        fnameStruct(1).fname_mani = 'C:\Neuroscience\imaging\FINAL\T7\E22_minLeaves_500_lmf_1_manifold_T7_FINAL.mat';
        %fnameStruct(1).fname_mani = 'C:\Neuroscience\imaging\E22_minLeaves_500_lmf_1_dim_2to7_manifold_T7_STEPSERROR.mat';
        
        fnameStruct(2).fname      = 'C:\Neuroscience\imaging\FINAL\E52_20171213_60per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(2).fname_mani = 'C:\Neuroscience\imaging\FINAL\T7\E52_minLeaves_500_lmf_1_manifold_T7_FINAL.mat';
        
        fnameStruct(3).fname      = 'C:\Neuroscience\imaging\FINAL\E53_20171214_60per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(3).fname_mani = 'C:\Neuroscience\imaging\FINAL\T7\E53_minLeaves_500_lmf_1_manifold_T7_FINAL.mat';
        
        fnameStruct(4).fname      = 'C:\Neuroscience\imaging\FINAL\E84_20190520_50per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(4).fname_mani = 'C:\Neuroscience\imaging\FINAL\T7\E84_minLeaves_500_lmf_1_manifold_T7_FINAL.mat';
        %fnameStruct(4).fname_mani = 'C:\Neuroscience\imaging\E84_minLeaves_500_lmf_1_dim_2to7_manifold_T7_20190520.mat';
        
    elseif strcmp(taskType, 'alternation')==1 || strcmp(taskType, 'Alternation')==1
        fnameStruct(1).fname      = 'C:\Neuroscience\imaging\FINAL\E38_20171020_50per_userSetSD11minDur0.modelingFINAL.mat';
        fnameStruct(1).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E38_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
        fnameStruct(2).fname      = 'C:\Neuroscience\imaging\FINAL\E44_20180118_50per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(2).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E44_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
        fnameStruct(3).fname      = 'C:\Neuroscience\imaging\FINAL\E48_20171024_60per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(3).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E48_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
        fnameStruct(4).fname      = 'C:\Neuroscience\imaging\FINAL\E62_20180214_70per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(4).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E62_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
        fnameStruct(5).fname      = 'C:\Neuroscience\imaging\FINAL\E63_20180323_50per_userSetSD11minDur0.modelingFINAL.mat';
        fnameStruct(5).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E63_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
        fnameStruct(6).fname      = 'C:\Neuroscience\imaging\FINAL\E66_20180208_50per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(6).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E66_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
        fnameStruct(7).fname      = 'C:\Neuroscience\imaging\FINAL\E67_20180326_50per_userSetSD11minDur0.modelingFINAL.mat';
        fnameStruct(7).fname_mani = 'C:\Neuroscience\imaging\FINAL\Alternation\E67_minLeaves_500_lmf_1_manifold_Alternation_FINAL.mat';
        
    elseif strcmp(taskType, 'alternationJeff')==1 || strcmp(taskType, 'AlternationJeff')==1
        fnameStruct(1).fname      = 'C:\Neuroscience\imaging\FINAL\E63_20180611_30per_userSetSD11minDur0.modelingFINAL.mat';
        fnameStruct(1).fname_mani = 'C:\Neuroscience\imaging\FINAL\AlternationJeff\E63_minLeaves_500_lmf_1_dim_2to7_manifold_AlternationJeff_FINAL.mat';
        
        fnameStruct(2).fname      = 'C:\Neuroscience\imaging\FINAL\E64_20180628_60per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(2).fname_mani = 'C:\Neuroscience\imaging\FINAL\AlternationJeff\E64_minLeaves_500_lmf_1_dim_2to7_manifold_AlternationJeff_FINAL.mat';
        
        fnameStruct(3).fname      = 'C:\Neuroscience\imaging\FINAL\E66_20180613_30per_userSetSD5minDur0.modelingFINAL.mat';
        fnameStruct(3).fname_mani = 'C:\Neuroscience\imaging\FINAL\AlternationJeff\E66_minLeaves_500_lmf_1_dim_2to7_manifold_AlternationJeff_FINAL.mat';
        
    end
    
elseif strcmp(user1, 'Manuel')==1
    
    if strcmp(taskTypes, 'towers')==1
        
        
        
    end
end
end