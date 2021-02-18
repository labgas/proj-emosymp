%% LaBGAS_copy_con_imgs_CANlab_folder_struct
%
% This simple script copies con2-4 images (pos, neu, neg) from subject/results folders in the root directory to
% subject folders in \fMRI_emotion_Giao\CANlab_second_level_analysis\data
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   November, 2020
%__________________________________________________________________________
% @(#)% LaBGAS_copy_con_imgs_CANlab_folder_struct     v1.1        
% last modified: 2021/02/18

%% define directories
rootdir='C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao';
canlabdatadir=fullfile(rootdir,'CANlab_second_level_analysis\data');

%% get subject folder names from rootdir
cd(rootdir);
hcs=dir('HC_*');
idh=[hcs.isdir]';
hcs={hcs(idh).name}';
pts=dir('PT_*');
idp=[pts.isdir]';
pts={pts(idp).name}';
subjs=[hcs;pts];

%% define anonymous function sm to write folder structure with cell array of folder names as input
sm=@(x)spm_mkdir(x); % defines spm_mkdir as an anonymous function sm

%% write folder structure in \fMRI_emotion_Giao\CANlab_second_level_analysis\data
cd(canlabdatadir);
cellfun(sm,subjs); % applies function sm to all cells of subjs

%% loop over subjects to copy the relevant con images over
for i=1:size(subjs,1)
    subjfirstleveldir=fullfile(rootdir,subjs{i},'results');
    subjcanlabdatadir=fullfile(canlabdatadir,subjs{i});
    cd(subjfirstleveldir);
    con_images=ls('con*.nii');
    con_images=con_images(2:4,:);
        for j=1:size(con_images,1)
            copyfile(fullfile(subjfirstleveldir,con_images(j,:)),fullfile(subjcanlabdatadir,con_images(j,:)));
        end
end