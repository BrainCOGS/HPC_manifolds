
%% Load the video and luminance data

load('C:\Neuroscience\imaging\FINAL\taskmanifold_Data\TT_embedding_rawdata.mat')
load('C:\Neuroscience\imaging\FINAL\taskmanifold_Data\video_luminance_rawdata.mat')


%% Run the code

fnameStruct = mind_makeFnameStruct('Edward','towers','laptop');

fitManifold2Video(fnameStruct(7).fname, '', './TT_embedding_rawdata.mat', '/Users/ms81/TT_embedding_manifold.mat', false)
