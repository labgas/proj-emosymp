%% LaBGAS_create_first_level_CANlab_tools_folder_structure
%
% This simple script 
% 1) creates the folder structure for first level analysis using CANlab tools
% 2) copies noise regressor and onset files generated by
% LaBGAS_extract_confound_reg_fMRIprep.m to run-specific directories in 
% 'derivatives\fmriprep' as this is consistent 
% with organization required by LaBGAS_get_firstlvl_dsgn_obj and
% LaBGAS_1_spm_fit_firstlvl_models
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   October, 2020
%__________________________________________________________________________
% @(#)% LaBGAS_create_first_level_CANlab_tools_folder_structure.m     v1.1       
% last modified: 2021/02/09

%% define directories (absolute path) and make new first level directory with model subdirectory
% comment out mkdir lines if you already manually created these dirs
rootdir='C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS';
% cd(rootdir);
% mkdir('firstlevel');
firstleveldir=fullfile(rootdir,'firstlevel');
% cd(firstleveldir);
% mkdir('model_1_CANlab_classic_GLM');
modeldir=fullfile(firstleveldir,'model_1_CANlab_classic_GLM');
derivdir=fullfile(rootdir,'derivatives\fmriprep');
runnames={'run-1';'run-2';'run-3';'run-4';'run-5';'run-6'}; % define names of run directories to be written in derivdir

%% get subject directory names from derivdir
cd(derivdir);
subjs=dir('sub-*');
idx=[subjs.isdir]';
subjs={subjs(idx).name}';

%% define anonymous function sm to write folder structure with cell array of directory names as input
sm=@(x)spm_mkdir(x); % defines spm_mkdir as an anonymous function sm

%% write folder structure in new model directory
cd(modeldir);
cellfun(sm,subjs); % applies function sm to all cells of subjs

%% loop over subjects to copy and move the onsets, noise_regs, functional scans, and fmriprep confound files to subdirectories per run
for i=54:size(subjs,1)
    subjmodeldir=fullfile(modeldir,subjs{i});
    subjderivdir=fullfile(derivdir,subjs{i},'func');
    cd(subjderivdir);
    onset_files=ls('onset*.mat');
    noise_files=ls('noise*.txt');
    nr_noise_files = size(noise_files,1);
    nr_onset_files = size(onset_files,1);
    cellfun(sm,runnames);
    niigz_files=cellstr(ls('s6*.nii.gz'));
    fmriprepconf_files=cellstr(ls('sub*.tsv'));
    if ~isequal(nr_noise_files,nr_onset_files,size(niigz_files,1),size(fmriprepconf_files,1))
        error('number of onset, noise, functional, and fmriprep confound files does not match');
    else
        for j=1:size(niigz_files,1)
            gunzip(niigz_files{j});
            niifile=ls('s6*.nii');
                if contains(niifile,runnames{j})
                    rundir=fullfile(subjderivdir,runnames{j});
                    movefile(fullfile(subjderivdir,niifile),fullfile(rundir,niifile));
                    movefile(fullfile(subjderivdir,onset_files(j,:)),fullfile(rundir,onset_files(j,:)));
                    movefile(fullfile(subjderivdir,noise_files(j,:)),fullfile(rundir,noise_files(j,:)));
                    copyfile(fullfile(subjderivdir,fmriprepconf_files{j,:}),fullfile(rundir,fmriprepconf_files{j,:}));
                    delete(niigz_files{j});
                elseif contains(niifile,runnames{j+1})
                    warning(strcat(runnames{j},'_missing_in_',subjs{i}));
                    rundir=fullfile(subjderivdir,runnames{j+1});
                    movefile(fullfile(subjderivdir,niifile),fullfile(rundir,niifile));
                    movefile(fullfile(subjderivdir,onset_files(j,:)),fullfile(rundir,onset_files(j,:)));
                    movefile(fullfile(subjderivdir,noise_files(j,:)),fullfile(rundir,noise_files(j,:)));
                    copyfile(fullfile(subjderivdir,fmriprepconf_files{j,:}),fullfile(rundir,fmriprepconf_files{j,:}));
                    delete(niigz_files{j});
                elseif contains(niifile,runnames{j+2})
                    warning(strcat(runnames{j+1},'_missing_in_',subjs{i})); % warnings may still be partly off - double-check whether missing runs match trouble in paradise
                    rundir=fullfile(subjderivdir,runnames{j+2});
                    movefile(fullfile(subjderivdir,niifile),fullfile(rundir,niifile));
                    movefile(fullfile(subjderivdir,onset_files(j,:)),fullfile(rundir,onset_files(j,:)));
                    movefile(fullfile(subjderivdir,noise_files(j,:)),fullfile(rundir,noise_files(j,:)));
                    copyfile(fullfile(subjderivdir,fmriprepconf_files{j,:}),fullfile(rundir,fmriprepconf_files{j,:}));
                    delete(niigz_files{j});
                elseif contains(niifile,runnames{j+3})
                    warning(strcat(runnames{j+2},'_missing_in_',subjs{i}));
                    rundir=fullfile(subjderivdir,runnames{j+3});
                    movefile(fullfile(subjderivdir,niifile),fullfile(rundir,niifile));
                    movefile(fullfile(subjderivdir,onset_files(j,:)),fullfile(rundir,onset_files(j,:)));
                    movefile(fullfile(subjderivdir,noise_files(j,:)),fullfile(rundir,noise_files(j,:)));
                    copyfile(fullfile(subjderivdir,fmriprepconf_files{j,:}),fullfile(rundir,fmriprepconf_files{j,:}));
                    delete(niigz_files{j});
                elseif contains(niifile,runnames{j+4})
                    warning(strcat(runnames{j+4},'_missing_in_',subjs{i}));
                    rundir=fullfile(subjderivdir,runnames{j+4});
                    movefile(fullfile(subjderivdir,niifile),fullfile(rundir,niifile));
                    movefile(fullfile(subjderivdir,onset_files(j,:)),fullfile(rundir,onset_files(j,:)));
                    movefile(fullfile(subjderivdir,noise_files(j,:)),fullfile(rundir,noise_files(j,:)));
                    copyfile(fullfile(subjderivdir,fmriprepconf_files{j,:}),fullfile(rundir,fmriprepconf_files{j,:}));
                    delete(niigz_files{j});
                end
        end
    end
end