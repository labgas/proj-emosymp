%% LaBGAS_copy_con_imgs_CANlab_folder_struct
%
% This simple script copies con1-3 images (pos, neu, neg) from
% rootdir\firstlevel\model_x_yyy\sub-zz\ to
% rootdir\secondlevel\model_x_yyy\con_images\sub-zz
% to prepare second level analysis with CANlab batch scripts
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   November, 2020
%__________________________________________________________________________
% @(#)% LaBGAS_copy_con_imgs_CANlab_folder_struct     v1.1        
% last modified: 2021/02/18

%% define model and directories
conimgs2copy=[1:3]; % numbers of con images to copy
model='model_1_CANlab_classic_GLM';
rootdir='C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS';
firstleveldir=fullfile(rootdir,'firstlevel\',model);
secondleveldir=fullfile(rootdir,'secondlevel\',model);
conimgsdir=fullfile(secondleveldir,'\data');

%% get subject directory names from firstleveldir
cd(firstleveldir);
subjs = dir([firstleveldir,'\sub-*']);
subjs = {subjs(:).name}';

%% define anonymous function sm to write folder structure with cell array of folder names as input
sm=@(x)spm_mkdir(x); % defines spm_mkdir as an anonymous function sm

%% write folder structure in conimgsdir
cd(conimgsdir);
cellfun(sm,subjs); % applies function sm to all cells of subjs

%% loop over subjects to copy the relevant con images over
for i=1:size(subjs,1)
    subjfirstleveldir=fullfile(firstleveldir,subjs{i});
    subjconimgsdir=fullfile(conimgsdir,subjs{i});
    cd(subjfirstleveldir);
    conimgs=ls('con*.nii');
    conimgs=conimgs(conimgs2copy,:);
        for j=1:size(conimgs,1)
            copyfile(fullfile(subjfirstleveldir,conimgs(j,:)),fullfile(subjconimgsdir,conimgs(j,:)));
        end
end