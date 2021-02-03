
%% Load the repository
addpath(genpath('C:\Edward\School\Princeton\PNI Code\HPC_manifolds\'));

%% Create structure of all behavioral logs for animals (at least 60% in T11 trials)

fname_E22_20170215 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170215.mat';
fname_E22_20170216 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170216.mat';
fname_E22_20170217 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170217.mat';
fname_E22_20170221 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170221.mat';
fname_E22_20170222 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170222.mat';
fname_E22_20170223 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170223.mat';
fname_E22_20170224 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170224.mat';
fname_E22_20170227 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170227.mat';
fname_E22_20170228 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170228.mat';
fname_E22_20170308 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170308.mat';
fname_E22_20170309 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170309.mat';
fname_E22_20170320 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170320.mat';
fname_E22_20170321 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170321.mat';
fname_E22_20170403 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170403.mat';
fname_E22_20170412 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_1\data\E22\PoissonBlocksReboot_cohort1_Bezos3_E22_T_20170412.mat';

outputBehavior_E22_20170215 = extractBehavior(fname_E22_20170215); 
outputBehavior_E22_20170216 = extractBehavior(fname_E22_20170216); 
outputBehavior_E22_20170217 = extractBehavior(fname_E22_20170217); 
outputBehavior_E22_20170221 = extractBehavior(fname_E22_20170221); 
outputBehavior_E22_20170222 = extractBehavior(fname_E22_20170222); 
outputBehavior_E22_20170223 = extractBehavior(fname_E22_20170223); 
outputBehavior_E22_20170224 = extractBehavior(fname_E22_20170224); 
outputBehavior_E22_20170227 = extractBehavior(fname_E22_20170227); 
outputBehavior_E22_20170228 = extractBehavior(fname_E22_20170228); 
outputBehavior_E22_20170308 = extractBehavior(fname_E22_20170308); 
outputBehavior_E22_20170309 = extractBehavior(fname_E22_20170309); 
outputBehavior_E22_20170320 = extractBehavior(fname_E22_20170320); 
outputBehavior_E22_20170321 = extractBehavior(fname_E22_20170321); 
outputBehavior_E22_20170403 = extractBehavior(fname_E22_20170403); 
outputBehavior_E22_20170412 = extractBehavior(fname_E22_20170412); 


fname_E39_20171011 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171011.mat';
fname_E39_20171012 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171012.mat';
fname_E39_20171102 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171102.mat';
fname_E39_20171103 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171103.mat';
fname_E39_20171110 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E39\PoissonBlocksReboot3_cohort3_Bezos3_E39_T_20171110.mat';

outputBehavior_E39_20171011 = extractBehavior(fname_E39_20171011);
outputBehavior_E39_20171012 = extractBehavior(fname_E39_20171012); 
outputBehavior_E39_20171102 = extractBehavior(fname_E39_20171102);
outputBehavior_E39_20171103 = extractBehavior(fname_E39_20171103);
outputBehavior_E39_20171110 = extractBehavior(fname_E39_20171110);


fname_E43_20170720 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170720.mat';
fname_E43_20170721 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170721.mat';
fname_E43_20170724 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170724.mat';
fname_E43_20170726 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170726.mat';
fname_E43_20170727 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170727.mat';
fname_E43_20170801 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170801.mat';
fname_E43_20170802 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170802.mat';
fname_E43_20170803 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170803.mat';
fname_E43_20170804 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170804.mat';
fname_E43_20170807 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170807.mat';
fname_E43_20170811 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170811.mat';
fname_E43_20170815 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170815.mat';
fname_E43_20170816 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170816.mat';
fname_E43_20170817 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170817.mat';
fname_E43_20170818 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170818.mat';
fname_E43_20170821 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170821.mat';
fname_E43_20170825 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170825.mat';
fname_E43_20170829 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170829.mat';
fname_E43_20170831 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170831.mat';
fname_E43_20170905 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170905.mat';
fname_E43_20170911 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170911.mat';
fname_E43_20170912 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E43\PoissonBlocksReboot2_cohort2_Bezos3_E43_T_20170912.mat';

outputBehavior_E43_20170720 = extractBehavior(fname_E43_20170720); 
outputBehavior_E43_20170721 = extractBehavior(fname_E43_20170721); 
outputBehavior_E43_20170724 = extractBehavior(fname_E43_20170724); 
outputBehavior_E43_20170726 = extractBehavior(fname_E43_20170726); 
outputBehavior_E43_20170727 = extractBehavior(fname_E43_20170727); 
outputBehavior_E43_20170801 = extractBehavior(fname_E43_20170801); 
outputBehavior_E43_20170802 = extractBehavior(fname_E43_20170802); 
outputBehavior_E43_20170803 = extractBehavior(fname_E43_20170803); 
outputBehavior_E43_20170804 = extractBehavior(fname_E43_20170804); 
outputBehavior_E43_20170807 = extractBehavior(fname_E43_20170807); 
outputBehavior_E43_20170811 = extractBehavior(fname_E43_20170811); 
outputBehavior_E43_20170815 = extractBehavior(fname_E43_20170815); 
outputBehavior_E43_20170816 = extractBehavior(fname_E43_20170816); 
outputBehavior_E43_20170817 = extractBehavior(fname_E43_20170817); 
outputBehavior_E43_20170818 = extractBehavior(fname_E43_20170818);
outputBehavior_E43_20170821 = extractBehavior(fname_E43_20170821);
outputBehavior_E43_20170825 = extractBehavior(fname_E43_20170825); 
outputBehavior_E43_20170829 = extractBehavior(fname_E43_20170829); 
outputBehavior_E43_20170831 = extractBehavior(fname_E43_20170831); 
outputBehavior_E43_20170905 = extractBehavior(fname_E43_20170905); 
outputBehavior_E43_20170911 = extractBehavior(fname_E43_20170911); 
outputBehavior_E43_20170912 = extractBehavior(fname_E43_20170912); 


fname_E44_20170914 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170914.mat';
fname_E44_20170915 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170915.mat';
fname_E44_20170918 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170918.mat';
fname_E44_20170919 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170919.mat';
fname_E44_20170920 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170920.mat';
fname_E44_20170921 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170921.mat';
fname_E44_20170922 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170922.mat';
fname_E44_20170925 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20170925.mat';
fname_E44_20171010 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171010.mat';
fname_E44_20171011 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171011.mat';
fname_E44_20171012 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171012.mat';
fname_E44_20171013 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171013.mat';
fname_E44_20171017 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171017.mat';
fname_E44_20171018 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171018.mat';
fname_E44_20171019 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171019.mat';
fname_E44_20171020 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171020.mat';
fname_E44_20171024 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171024.mat';
fname_E44_20171025 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171025.mat';
fname_E44_20171108 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_2\data\E44\PoissonBlocksReboot2_cohort2_Bezos3_E44_T_20171108.mat';

outputBehavior_E44_20170914 = extractBehavior(fname_E44_20170914); 
outputBehavior_E44_20170915 = extractBehavior(fname_E44_20170915); 
outputBehavior_E44_20170918 = extractBehavior(fname_E44_20170918); 
outputBehavior_E44_20170919 = extractBehavior(fname_E44_20170919); 
outputBehavior_E44_20170920 = extractBehavior(fname_E44_20170920);
outputBehavior_E44_20170921 = extractBehavior(fname_E44_20170921); 
outputBehavior_E44_20170922 = extractBehavior(fname_E44_20170922); 
outputBehavior_E44_20170925 = extractBehavior(fname_E44_20170925);
outputBehavior_E44_20171010 = extractBehavior(fname_E44_20171010); 
outputBehavior_E44_20171011 = extractBehavior(fname_E44_20171011);
outputBehavior_E44_20171012 = extractBehavior(fname_E44_20171012); 
outputBehavior_E44_20171013 = extractBehavior(fname_E44_20171013); 
outputBehavior_E44_20171017 = extractBehavior(fname_E44_20171017); 
outputBehavior_E44_20171018 = extractBehavior(fname_E44_20171018); 
outputBehavior_E44_20171019 = extractBehavior(fname_E44_20171019); 
outputBehavior_E44_20171020 = extractBehavior(fname_E44_20171020); 
outputBehavior_E44_20171024 = extractBehavior(fname_E44_20171024); 
outputBehavior_E44_20171025 = extractBehavior(fname_E44_20171025); 
outputBehavior_E44_20171108 = extractBehavior(fname_E44_20171108); 


fname_E47_20170927 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20170927.mat';
fname_E47_20171003 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171003.mat';
fname_E47_20171004 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171004.mat';
fname_E47_20171005 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171005.mat';
fname_E47_20171006 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171006.mat';
fname_E47_20171010 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171010.mat';
fname_E47_20171011 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171011.mat';
fname_E47_20171012 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171012.mat';
fname_E47_20171013 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171013.mat';
fname_E47_20171018 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E47\PoissonBlocksReboot3_cohort3_Bezos3_E47_T_20171018.mat';

outputBehavior_E47_20170927 = extractBehavior(fname_E47_20170927); 
outputBehavior_E47_20171003 = extractBehavior(fname_E47_20171003); 
outputBehavior_E47_20171004 = extractBehavior(fname_E47_20171004); 
outputBehavior_E47_20171005 = extractBehavior(fname_E47_20171005); 
outputBehavior_E47_20171006 = extractBehavior(fname_E47_20171006); 
outputBehavior_E47_20171010 = extractBehavior(fname_E47_20171010); 
outputBehavior_E47_20171011 = extractBehavior(fname_E47_20171011);
outputBehavior_E47_20171012 = extractBehavior(fname_E47_20171012); 
outputBehavior_E47_20171013 = extractBehavior(fname_E47_20171013);
outputBehavior_E47_20171018 = extractBehavior(fname_E47_20171018); 


fname_E48_20170807 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170807.mat';
fname_E48_20170808 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170808.mat';
fname_E48_20170809 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170809.mat';
fname_E48_20170810 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170810.mat';
fname_E48_20170811 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170811.mat';
fname_E48_20170814 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170814.mat';
fname_E48_20170815 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170815.mat';
fname_E48_20170816 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170816.mat';
fname_E48_20170817 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170817.mat';
fname_E48_20170823 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170823.mat';
fname_E48_20170824 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170824.mat';
fname_E48_20170828 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170828.mat';
fname_E48_20170829 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170829.mat';
fname_E48_20170830 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170830.mat';
fname_E48_20170831 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170831.mat';
fname_E48_20170905 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170905.mat';
fname_E48_20170906 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170906.mat';
fname_E48_20170907 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170907.mat';
fname_E48_20170908 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170908.mat';
fname_E48_20170911 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170911.mat';
fname_E48_20170914 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170914.mat';
fname_E48_20170915 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170915.mat';
fname_E48_20170918 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170918.mat';
fname_E48_20170919 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170919.mat';
fname_E48_20170920 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20170920.mat';
fname_E48_20171004 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_3\data\E48\PoissonBlocksReboot3_cohort3_Bezos3_E48_T_20171004.mat';

outputBehavior_E48_20170807 = extractBehavior(fname_E48_20170807); 
outputBehavior_E48_20170808 = extractBehavior(fname_E48_20170808);
outputBehavior_E48_20170809 = extractBehavior(fname_E48_20170809);
outputBehavior_E48_20170810 = extractBehavior(fname_E48_20170810);
outputBehavior_E48_20170811 = extractBehavior(fname_E48_20170811);
outputBehavior_E48_20170814 = extractBehavior(fname_E48_20170814); 
outputBehavior_E48_20170815 = extractBehavior(fname_E48_20170815);
outputBehavior_E48_20170816 = extractBehavior(fname_E48_20170816);
outputBehavior_E48_20170817 = extractBehavior(fname_E48_20170817);
outputBehavior_E48_20170823 = extractBehavior(fname_E48_20170823); 
outputBehavior_E48_20170824 = extractBehavior(fname_E48_20170824); 
outputBehavior_E48_20170828 = extractBehavior(fname_E48_20170828); 
outputBehavior_E48_20170829 = extractBehavior(fname_E48_20170829); 
outputBehavior_E48_20170830 = extractBehavior(fname_E48_20170830); 
outputBehavior_E48_20170831 = extractBehavior(fname_E48_20170831); 
outputBehavior_E48_20170905 = extractBehavior(fname_E48_20170905); 
outputBehavior_E48_20170906 = extractBehavior(fname_E48_20170906); 
outputBehavior_E48_20170907 = extractBehavior(fname_E48_20170907); 
outputBehavior_E48_20170908 = extractBehavior(fname_E48_20170908);
outputBehavior_E48_20170911 = extractBehavior(fname_E48_20170911); 
outputBehavior_E48_20170914 = extractBehavior(fname_E48_20170914); 
outputBehavior_E48_20170915 = extractBehavior(fname_E48_20170915); 
outputBehavior_E48_20170918 = extractBehavior(fname_E48_20170918); 
outputBehavior_E48_20170919 = extractBehavior(fname_E48_20170919); 
outputBehavior_E48_20170920 = extractBehavior(fname_E48_20170920); 
outputBehavior_E48_20171004 = extractBehavior(fname_E48_20171004); 


fname_E65_20180131 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180131.mat';
fname_E65_20180201 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180201.mat';
fname_E65_20180202 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180202.mat';
fname_E65_20180205 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180205.mat';
fname_E65_20180206 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180206.mat';
fname_E65_20180207 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180207.mat';
fname_E65_20180208 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180208.mat';
fname_E65_20180209 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180209.mat';
fname_E65_20180212 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180212.mat';
fname_E65_20180213 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180213.mat';
fname_E65_20180214 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180214.mat';
fname_E65_20180312 = '\\bucket.pni.princeton.edu\braininit\RigData\scope\bay3\edward\PoissonTowers_4\data\E65\PoissonBlocksReboot4_cohort4_Bezos3_E65_T_20180312.mat';

outputBehavior_E65_20180131 = extractBehavior(fname_E65_20180131); 
outputBehavior_E65_20180201 = extractBehavior(fname_E65_20180201); 
outputBehavior_E65_20180202 = extractBehavior(fname_E65_20180202); 
outputBehavior_E65_20180205 = extractBehavior(fname_E65_20180205); 
outputBehavior_E65_20180206 = extractBehavior(fname_E65_20180206); 
outputBehavior_E65_20180207 = extractBehavior(fname_E65_20180207); 
outputBehavior_E65_20180208 = extractBehavior(fname_E65_20180208);
outputBehavior_E65_20180209 = extractBehavior(fname_E65_20180209); 
outputBehavior_E65_20180212 = extractBehavior(fname_E65_20180212); 
outputBehavior_E65_20180213 = extractBehavior(fname_E65_20180213);
outputBehavior_E65_20180214 = extractBehavior(fname_E65_20180214); 
outputBehavior_E65_20180312 = extractBehavior(fname_E65_20180312); 


S = CatStructFields(2, outputBehavior_E22_20170215, ...
    outputBehavior_E22_20170216, ...
    outputBehavior_E22_20170217, ...
    outputBehavior_E22_20170221, ...
    outputBehavior_E22_20170222, ...
    outputBehavior_E22_20170223, ...
    outputBehavior_E22_20170224, ...
    outputBehavior_E22_20170227, ...
    outputBehavior_E22_20170228, ...
    outputBehavior_E22_20170308, ...
    outputBehavior_E22_20170309, ...
    outputBehavior_E22_20170320, ...
    outputBehavior_E22_20170321, ...
    outputBehavior_E22_20170403, ...
    outputBehavior_E22_20170412, ...
    outputBehavior_E39_20171011, ...
    outputBehavior_E39_20171012, ...
    outputBehavior_E39_20171102, ...
    outputBehavior_E39_20171103, ...
    outputBehavior_E39_20171110, ...
    outputBehavior_E43_20170720, ...
    outputBehavior_E43_20170721, ...
    outputBehavior_E43_20170724, ...
    outputBehavior_E43_20170726, ...
    outputBehavior_E43_20170727, ... 
    outputBehavior_E43_20170801, ...
    outputBehavior_E43_20170802, ...
    outputBehavior_E43_20170803, ...
    outputBehavior_E43_20170804, ...
    outputBehavior_E43_20170807, ...
    outputBehavior_E43_20170811, ...
    outputBehavior_E43_20170815, ...
    outputBehavior_E43_20170816, ...
    outputBehavior_E43_20170817, ...
    outputBehavior_E43_20170818, ...
    outputBehavior_E43_20170821, ...
    outputBehavior_E43_20170825, ...
    outputBehavior_E43_20170829, ...
    outputBehavior_E43_20170831, ...
    outputBehavior_E43_20170905, ...
    outputBehavior_E43_20170911, ...
    outputBehavior_E43_20170912, ...
    outputBehavior_E44_20170914, ...
    outputBehavior_E44_20170915, ...
    outputBehavior_E44_20170918, ...
    outputBehavior_E44_20170919, ...
    outputBehavior_E44_20170920, ...
    outputBehavior_E44_20170921, ...
    outputBehavior_E44_20170922, ...
    outputBehavior_E44_20170925, ...
    outputBehavior_E44_20171010, ...
    outputBehavior_E44_20171011, ...
    outputBehavior_E44_20171012, ...
    outputBehavior_E44_20171013, ...
    outputBehavior_E44_20171017, ...
    outputBehavior_E44_20171018, ...
    outputBehavior_E44_20171019, ...
    outputBehavior_E44_20171020, ...
    outputBehavior_E44_20171024, ...
    outputBehavior_E44_20171025, ...
    outputBehavior_E44_20171108, ...
    outputBehavior_E47_20170927, ...
    outputBehavior_E47_20171003, ...
    outputBehavior_E47_20171004, ...
    outputBehavior_E47_20171005, ...
    outputBehavior_E47_20171006, ...
    outputBehavior_E47_20171010, ...
    outputBehavior_E47_20171011, ...
    outputBehavior_E47_20171012, ...
    outputBehavior_E47_20171013, ...
    outputBehavior_E47_20171018, ...
    outputBehavior_E48_20170807, ...
    outputBehavior_E48_20170808, ...
    outputBehavior_E48_20170809, ...
    outputBehavior_E48_20170810, ...
    outputBehavior_E48_20170811, ...
    outputBehavior_E48_20170814, ...
    outputBehavior_E48_20170815, ...
    outputBehavior_E48_20170816, ...
    outputBehavior_E48_20170817, ...
    outputBehavior_E48_20170823, ...
    outputBehavior_E48_20170824, ...
    outputBehavior_E48_20170828, ...
    outputBehavior_E48_20170829, ... 
    outputBehavior_E48_20170830, ...
    outputBehavior_E48_20170831, ...
    outputBehavior_E48_20170905, ...
    outputBehavior_E48_20170906, ...
    outputBehavior_E48_20170907, ...
    outputBehavior_E48_20170908, ...
    outputBehavior_E48_20170911, ...
    outputBehavior_E48_20170914, ...
    outputBehavior_E48_20170915, ...
    outputBehavior_E48_20170918, ...
    outputBehavior_E48_20170919, ...
    outputBehavior_E48_20170920, ...
    outputBehavior_E48_20171004, ...
    outputBehavior_E65_20180131, ...
    outputBehavior_E65_20180201, ...
    outputBehavior_E65_20180202, ...
    outputBehavior_E65_20180205, ...
    outputBehavior_E65_20180206, ...
    outputBehavior_E65_20180207, ...
    outputBehavior_E65_20180208, ...
    outputBehavior_E65_20180209, ...
    outputBehavior_E65_20180212, ...
    outputBehavior_E65_20180213, ...
    outputBehavior_E65_20180214, ...
    outputBehavior_E65_20180312);


%% Make the logistic regression structure

logisticRegr = logRegressionFromConcatLog2(S);

%% Make the figures using code copied from Lucas's owf_taskPerfComp

layout       = [ 1 2 ];
fig          = PaneledFigure(layout, 'smaller');
iPanel       = 0;

%% psychometrics for all mice on the paper
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs, 'on')
for iMouse = 1:logisticRegr.nGoodMice
    plot(axs, toBinCenters(logisticRegr.bins),                    ...
        logisticRegr.goodValues(:,iMouse),                  ...
        '-','linewidth',  FormatDefaults.linewidthThin,     ...
        'color',      FormatDefaults.lightGray);
end
errorbar(axs, toBinCenters(logisticRegr.bins), logisticRegr.values, logisticRegr.sem, ...
    'k-','linewidth',FormatDefaults.linewidthRegular, 'markersize', 4);
xlim([0 200]); ylim([-.1 .4])
set(axs, 'xtick', [0 100 200], 'ytick', 0:.1:.4)
xlabel('Cue y (cm)')
ylabel('Weight on decision')
psychByMouse = xMousePyschometric(S.choice,S.nCues_RminusL,S.mouseID);
psych        = psychometricFit(S.choice,S.nCues_RminusL,true);

%% plot best fitting sigmoid for each mouse and overlay metamouse
iPanel                  = iPanel + 1;
axs                     = fig.panel(iPanel);
hold(axs, 'on')
for iMouse = 1:size(psychByMouse.goodCurves,2)
    plot(psychByMouse.fitXaxis,100.*psychByMouse.goodFits(:,iMouse),'-', ...
        'linewidth',         FormatDefaults.linewidthThin,     ...
        'color'    ,         FormatDefaults.lightGray)
end
plotPsychometricCurve_ctrl(psych, 'all', axs, 'k', true, 'k', false);
set(axs, 'ytick', 0:25:100)
ylabel(axs,'Went right (%)')
xlabel(axs, '\Delta towers (#R - #L)')

%% Length of Trial Durations

outputTrialDurations_E22 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E22_20170227_30per_userSetSD11minDur0.modelingFINAL.mat');
outputTrialDurations_E39 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E39_20171103_40per_userSetSD11minDur0.modelingFINAL.mat');
outputTrialDurations_E43 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E43_20170802_70per_userSetSD5minDur0.modelingFINAL.mat');
outputTrialDurations_E44 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E44_20171018_50per_userSetSD5minDur0.modelingFINAL.mat');
outputTrialDurations_E47 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E47_20170927_70per_userSetSD5minDur0.modelingFINAL.mat');
outputTrialDurations_E48 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E48_20170829_70per_userSetSD5minDur0.modelingFINAL.mat');
outputTrialDurations_E65 = extractTrialDurations('C:\Neuroscience\imaging\FINAL\E65_20180202_60per_userSetSD5minDur0.modelingFINAL.mat');





