%% LaBGAS_create_single_trial_fmri_data_st_obj
%
% This script creates an fmri_data_st object from the single trial con
% images written in rootdir\secondlevel\model_x_yyy\data\sub-zz by
% LaBGAS_copy_single_trial_con_imgs_CANlab_folder_struct.m and adds a
% convenient metadata_table field containing
% 1. single trial ratings from rootdir\phenotype\phenotype_trial.tsv
% 2. single trial vifs also copied by the abovementioned script
%
% Con images include con_0007:end defined and created by 
% LaBGAS_add_single_trial_contrasts.m
% These are 2 con images per condition per run, for each single trial,
% ordered per condition: negative, neutral, positive
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   March, 2021
%__________________________________________________________________________
% @(#)% LaBGAS_create_single_trial_fmri_data_st_obj     v1.0        
% last modified: 2021/03/30


%% define directories

model = 'model_2_CANlab_single_trial_GLM';
rootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS';
phenodir = fullfile(rootdir,'phenotype');
secondleveldir = fullfile(rootdir,'secondlevel\',model);
resultsdir = fullfile(secondleveldir,'results');
conimgsdir = fullfile(secondleveldir,'\data');


%% get subject directory names from secondleveldir, and behavioral data from phenodir

load(fullfile(resultsdir,'image_names_and_setup.mat'));
subjs = ls([conimgsdir,'\sub-*']);
idx_behav1 = ~isnan(DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu); % index for subjects with missing behavioral data, which we want to exclude
subjs = subjs(idx_behav1,:);


%% read con images and vifs, and add to fmri_data_st object

conimgs = [];
vif_values = [];
vif_names = [];
for i = 1:size(subjs,1)
    subjdir = fullfile(conimgsdir,subjs(i,:));
    subjconimgs = ls([subjdir,'\con_*']);
    subjconimgs = subjconimgs(4:end,:);
    for j = 1:size(subjconimgs,1)
        subjconimgs_fp{j,1} = [deblank(subjdir),'\',subjconimgs(j,:)];
    end
    conimgs = [conimgs;subjconimgs_fp];
    clear subjconimgs subjconimgs_fp
    load(fullfile(deblank(subjdir),'vifs.mat'));
    vif_values = [vif_values;vifs.allvifs(1:end-1)'];
    vif_names = [vif_names;char(vifs.name(1:end-1)')];
    clear vifs
end

this_dat = fmri_data_st(conimgs);
this_dat = remove_empty(this_dat);
this_dat.source_notes = 'single trial con images for emosymp dataset';
this_dat.mask_descrip = 'default CANlab brainmask';
this_dat.metadata_table.vif_values = vif_values;
this_dat.metadata_table.vif_names = vif_names;


%% add ratings to metadata_table field of fmri_data_st object

behav_dat = readtable(fullfile(phenodir,'phenotype_trial.tsv'),'Filetype','text','Delimiter','tab');
behav_dat.patient(behav_dat.patient == 0) = -1; % reference coding (-1 1) rather than effects coding (0 1) for reasons of consistency with other scripts
idx_participant = logical(behav_dat.participant_included);
idx_run = logical(behav_dat.run_included);
idx_behav2 = ~isnan(behav_dat.valencepresentation_within_run);
idx = idx_participant & idx_run & idx_behav2;
behav_dat = behav_dat(idx,:);
behav_dat = sortrows(behav_dat,{'subject_id','trial_valence'});
this_dat.metadata_table = [this_dat.metadata_table behav_dat];


%% add behavioral outcome you want to predict to Y field of fmri_data_st object

this_dat.Y = this_dat.metadata_table.state_symptoms;
this_dat.Y_descrip = 'somatic symptom rating';

%% save fmri_data_st object

savefilename = fullfile(resultsdir, 'single_trial_fmri_data_st_object.mat');
save(savefilename, 'this_dat');