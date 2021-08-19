% This script runs single-level mediation models for proj-emosymp 
% (Van den Houte, Bogaerts, et al, in prep), with the following variables
% X: group (patients versus controls)
% M: different options
% 1. NPS response (including subregions) for each contrast
% 2. voxel-wise brain mediator
% 3. Principal Directions of Mediation (PDM) multivariate brain mediator
% Y: somatic symptom ratings for each contrast
%
% It requires all study-specific and most core prep_ scripts for proj-emosymp 
% (https://github.com/labgas/proj-emosymp) and CANlab_help examples (LaBGAS
% fork, https://github.com/labgas/CANlab_help_examples) 
% to be run first, or load your saved .mat files if you ran these
% scripts earlier!
%
% Options for masking and scaling are currently still hard-coded, but 
% can easily be set as options in a2_emosymp_m1_s2_set_default_options, 
% or as part of this script at a later stage if desired
%
% This script is adapted from the following CANlab scripts:
%
% https://www.dropbox.com/s/uvfhdstfveve393/s8_mediation_model_FAPS.m?dl=0 
% by @pkragel
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
% by @tor and @martin
% 
% see also 
% https://canlab.github.io/_pages/mediation_example_script_1/mediation_example_script_1.html
% https://canlab.github.io/_pages/mediation_brain_sample_report/mediation_brain_results_report.html
%
% @lukasvo76 @Dartmouth, March-April 2021
% 
% WORK IN PROGRESS
% FIGURES, TITLES, AND STRUCTURE OF HTML OUTPUT CAN BE IMPROVED
% LOOPING OVER SUBPATTERNS/REGIONS CAN BE IMPLEMENTED


%% DEFINE PATHS AND VARIABLES TO BE USED IN ALL ANALYSES
% -------------------------------------------------------------------------
addpath(genpath('C:\Users\lukas\Documents\GitHub\MediationToolbox')); % add CANlab mediation toolbox to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\proj-emosymp')); % add proj-emosymp scripts to your path

a_emosymp_m1_s1_set_up_paths_always_run_first

load(fullfile(resultsdir,'image_names_and_setup.mat'));

mediationdir = fullfile(resultsdir,'mediation_analysis'); % a_ script needs to be run first (always!) so resultsdir is known
brain_mediationdir = fullfile(mediationdir,'brain_voxelwise');
pdm_mediationdir = fullfile(mediationdir,'brain_pdm');
if ~isfolder(mediationdir)
    mkdir(mediationdir);
end
if ~isfolder(brain_mediationdir)
    mkdir(brain_mediationdir);
end
if ~isfolder(pdm_mediationdir)
    mkdir(pdm_mediationdir);
end
cd(mediationdir);

idx = ~isnan(DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu); % index for subjects with missing behavioral data, which we want to exclude
X = DAT.BETWEENPERSON.group; % define predictor
Ydat = [DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu, DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_pos, DAT.BEHAVIOR.behavioral_data_table.symptoms_pos_neu]; % define matrix of outcomes (subjects * contrasts)
Ynames = {'somatic_symptoms'};
contrastnames = DAT.SIG_contrasts.raw.dotproduct.conditionnames; % get names of contrasts


%% NPS
% -------------------------------------------------------------------------

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(:,i));
    Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(:,i));
    Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(:,i));
    
    % run mediation models for each of the mediators
    [paths, toplevelstats, ~] = mediation(X,Y{i},M{i},'names',{'group',contrastnames{i},'NPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPS{i} = paths;
    toplevelstats_NPS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mpos{i},'names',{'group',contrastnames{i},'NPSpos response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPSpos{i} = paths;
    toplevelstats_NPSpos{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mneg{i},'names',{'group',contrastnames{i},'NPSneg response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPSneg{i} = paths;
    toplevelstats_NPSneg{i} = toplevelstats;
    clear paths toplevelstats;
    
end


%% SELECTED NPS SUBREGIONS
% -------------------------------------------------------------------------

for i = 1:size(contrastnames,2) 
    
    Y{i} = Ydat(:,i);
    
    MrIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(:,2);
    MrdpIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(:,6);
    MrS2_Op{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(:,7);
    MdACC{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(:,8);
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MrIns{i},'names',{'group',contrastnames{i},'rIns response'},'boottop','plots');
    drawnow;snapnow;
    paths_rIns{i} = paths;
    toplevelstats_rIns{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MrdpIns{i},'names',{'group',contrastnames{i},'rdpIns response'},'boottop','plots');
    drawnow;snapnow;
    paths_rdpIns{i} = paths;
    toplevelstats_rdpIns{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MrS2_Op{i},'names',{'group',contrastnames{i},'rS2_Op response'},'boottop','plots');
    drawnow;snapnow;
    paths_rS2_Op{i} = paths;
    toplevelstats_rS2_Op{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MdACC{i},'names',{'group',contrastnames{i},'dACC response'},'boottop','plots');
    drawnow;snapnow;
    paths_dACC{i} = paths;
    toplevelstats_dACC{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},MrIns{i},'M',[MrdpIns{i},MrS2_Op{i},MdACC{i}],'names',{'group',contrastnames{i},'rIns response','rdpIns response','rS2_Op response','dACC response'},'boottop','plots');
    drawnow;snapnow;
    paths_mult_med{i} = paths;
    toplevelstats_mult_med{i} = toplevelstats;
    clear paths toplevelstats;
    
end


%% VOXEL-BASED MEDIATION
%--------------------------------------------------------------------------

for i = 1:size(contrastnames,2)
    
    brain_mediationdir_contrast = fullfile(brain_mediationdir, contrastnames{i});
    if ~isfolder(brain_mediationdir_contrast)
        mkdir(brain_mediationdir_contrast);
    end
    cd(brain_mediationdir_contrast);
    
    Y{i} = Ydat(:,i);
    
    M_brain{i} = DAT.imgs{i}(:);
    
    mediation_brain(X,Y{i},char(M_brain{i}),'mask',maskname_glm,'names',{'group',contrastnames{i},'brain'});
    
    cd(brain_mediationdir);
    
end

%% now run the following script from each results directory

publish_mediation_report;


%% MULTIVARIATE MEDIATION
%--------------------------------------------------------------------------
% based on
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
% and personal communication with Martin Lindquist and Xiaochun Han

X_c = num2cell(X(idx)); % convert double to cell array

for i = 1:size(contrastnames,2)
    
    pdm_mediationdir_contrast = fullfile(pdm_mediationdir, contrastnames{i});
    if ~isfolder(pdm_mediationdir_contrast)
        mkdir(pdm_mediationdir_contrast);
    end
    cd(pdm_mediationdir_contrast);
    
    pdm_mediationdir_contrast_behav = fullfile(pdm_mediationdir_contrast, Ynames{1});
            if ~isfolder(pdm_mediationdir_contrast_behav)
                mkdir(pdm_mediationdir_contrast_behav);
            end
    cd(pdm_mediationdir_contrast_behav);
    
    if ~isfile(strcat('PDMresults_',contrastnames{i},'_',Ynames{1},'.mat')) || ~isfile(strcat('data_objects_',contrastnames{i},'_',Ynames{1},'.mat'))
        
        Y{i} = Ydat(idx,i);
        Y_c{i} = num2cell(Y{i});

        M_brain{i} = DAT.imgs{i}(idx);
        names = M_brain{i};

        mask = which('gray_matter_mask.nii');

        for k=1:size(X_c,1)
            dat = fmri_data(names{k},mask); 
            m{k} = dat.dat; 
        end
        
        save(strcat('data_objects_',contrastnames{i},'_',Ynames{1},'.mat'),'dat','-v7.3');

        pdm = multivariateMediation(X_c,Y_c{i},m,'B',20,'svd','plots');
        pdmfull = pdm;
        pdm = multivariateMediation(pdm,'nPDM',2);
        pdm = multivariateMediation(pdm,'noPDMestimation','bootPDM',1:2,'bootjPDM','Bsamp',5000,'save2file',strcat('PDMresults_',contrastnames{i},'_',Ynames{1},'.mat'));
        save(strcat('PDMresults_',contrastnames{i},'_',Ynames{1},'.mat'),'pdmfull','-append');
    
    else
        
        load(strcat('PDMresults_',contrastnames{i},'_',Ynames{1},'.mat'));
        pdm = out;
        clear out;
        load(strcat('data_objects_',contrastnames{i},'_',Ynames{1},'.mat'));
        
    end % if loop
    
    % source reconstruction
        % "source reconstruction is thus no more than a covariance map
        % which shows how much each voxel covaries with model predictions
        % across images" - Bogdan Petre
        % ref Haufe et al NeuroImage 2014
        
    for n = 1:size(pdm.Wfull,2)

        X_source{i} = fmri_data(DAT.imgs{i}(idx),which('gray_matter_mask.nii'));
        Xz_source{i} = rescale(X_source{i},'centervoxels');
        M_source{i,n} = pdm.Wfull{1,n};
        P_source{i,n} = M_source{i,n}'*Xz_source{i}.dat;
        source{i,n} = (P_source{i,n}*Xz_source{i}.dat');
        source{i,n} = source{i,n}'./(size(P_source{i,n},2)-1);
        source_obj{i,n} = dat;
        source_obj{i,n}.dat = source{i,n};
        source_obj{i,n}.image_names = 'source reconstruction';
        source_obj{i,n}.fullpath = '';
        source_obj{i,n}.history = {['source reconstructed from con images and pdm_',num2str(n)]};
        source_obj_j{n} = source_obj{i,n};

    end % source recon loop

    save(strcat('PDM_source_recon_',contrastnames{i},'_',Ynames{1},'.mat'),'source_obj_j');
        
    cd(pdm_mediationdir);
    
end % for loop contrasts

%% now run the following script from each results directory

publish_multivariate_mediation_report;
