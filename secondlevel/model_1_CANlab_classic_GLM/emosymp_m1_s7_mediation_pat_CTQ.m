% This script runs single-level mediation models for proj-emosymp 
% (Van den Houte, Bogaerts, et al, in prep), in patients only, 
% with the following variables
% X: CTQ
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
% @lukasvo76 @KULeuven, May 2021
% 
% WORK IN PROGRESS
% FIGURES, TITLES, AND STRUCTURE OF HTML OUTPUT CAN BE IMPROVED
% LOOPING OVER SUBPATTERNS/REGIONS CAN BE IMPLEMENTED


%% DEFINE PATHS AND VARIABLES TO BE USED IN ALL ANALYSES
% -------------------------------------------------------------------------
addpath(genpath('C:\Users\lukas\Documents\GitHub\MediationToolbox')); % add CANlab mediation toolbox to your path

a_emosymp_m1_s1_set_up_paths_always_run_first

load(fullfile(resultsdir,'image_names_and_setup.mat'));

mediationdir = fullfile(resultsdir,'mediation_analysis'); % a_ script needs to be run first (always!) so resultsdir is known
ctq_brain_mediationdir = fullfile(mediationdir,'pat_ctq_brain_voxelwise');
ctq_pdm_mediationdir = fullfile(mediationdir,'pat_ctq_brain_pdm');
if ~isfolder(mediationdir)
    mkdir(mediationdir);
end
if ~isfolder(ctq_brain_mediationdir)
    mkdir(ctq_brain_mediationdir);
end
if ~isfolder(ctq_pdm_mediationdir)
    mkdir(ctq_pdm_mediationdir);
end
cd(mediationdir);

idx_ctq = ~isnan(DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu) & ~isnan(DAT.BEHAVIOR.behavioral_data_table.CTQ) & (DAT.BEHAVIOR.behavioral_data_table.patient == 1); % index for subjects with missing behavioral data, ctq data, and healthy controls, all of which we want to exclude
X = zscore(DAT.BEHAVIOR.behavioral_data_table.CTQ(idx_ctq)); % define predictor and zscore it
Ydat = [DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu(idx_ctq), DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_pos(idx_ctq), DAT.BEHAVIOR.behavioral_data_table.symptoms_pos_neu(idx_ctq)]; % define matrix of outcomes (subjects * contrasts)
contrastnames = DAT.SIG_contrasts.raw.dotproduct.conditionnames; % get names of contrasts


%% NPS
% -------------------------------------------------------------------------

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(idx_ctq,i));
    Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(idx_ctq,i));
    Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(idx_ctq,i));
    
    % run mediation models for each of the mediators
    [paths, toplevelstats, ~] = mediation(X,Y{i},M{i},'names',{'CTQ',contrastnames{i},'NPS response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPS{i} = paths;
    toplevelstats_NPS{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mpos{i},'names',{'CTQ',contrastnames{i},'NPSpos response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPSpos{i} = paths;
    toplevelstats_NPSpos{i} = toplevelstats;
    clear paths toplevelstats;
    
    [paths, toplevelstats, ~] = mediation(X,Y{i},Mneg{i},'names',{'CTQ',contrastnames{i},'NPSneg response'},'boottop','plots');
    drawnow;snapnow;
    paths_NPSneg{i} = paths;
    toplevelstats_NPSneg{i} = toplevelstats;
    clear paths toplevelstats;
    
end


% %% SELECTED NPS SUBREGIONS
% % -------------------------------------------------------------------------
% 
% for i = 1:size(contrastnames,2) 
%     
%     Y{i} = Ydat(:,i);
%     
%     MrIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx_ctq,2);
%     MrdpIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx_ctq,6);
%     MrS2_Op{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx_ctq,7);
%     MdACC{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx_ctq,8);
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrIns{i},'names',{'group',contrastnames{i},'rIns response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_rIns{i} = paths;
%     toplevelstats_rIns{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrdpIns{i},'names',{'group',contrastnames{i},'rdpIns response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_rdpIns{i} = paths;
%     toplevelstats_rdpIns{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrS2_Op{i},'names',{'group',contrastnames{i},'rS2_Op response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_rS2_Op{i} = paths;
%     toplevelstats_rS2_Op{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MdACC{i},'names',{'group',contrastnames{i},'dACC response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_dACC{i} = paths;
%     toplevelstats_dACC{i} = toplevelstats;
%     clear paths toplevelstats;
%     
%     [paths, toplevelstats, ~] = mediation(X,Y{i},MrIns{i},'M',[MrdpIns{i},MrS2_Op{i},MdACC{i}],'names',{'group',contrastnames{i},'rIns response','rdpIns response','rS2_Op response','dACC response'},'boottop','plots');
%     drawnow;snapnow;
%     paths_mult_med{i} = paths;
%     toplevelstats_mult_med{i} = toplevelstats;
%     clear paths toplevelstats;
%     
% end


%% VOXEL-BASED MEDIATION
%--------------------------------------------------------------------------

for i = 1:size(contrastnames,2)
    
    brain_mediationdir_contrast = fullfile(ctq_brain_mediationdir, contrastnames{i});
    if ~isfolder(brain_mediationdir_contrast)
        mkdir(brain_mediationdir_contrast);
    end
    cd(brain_mediationdir_contrast);
    
    Y{i} = Ydat(:,i);
    
    M_brain{i} = DAT.imgs{i}(idx_ctq);
    
    mediation_brain(X,Y{i},char(M_brain{i}),'mask',maskname_glm,'names',{'group',contrastnames{i},'brain'});
    
    cd(ctq_brain_mediationdir);
    
end

% now run the following script from each results directory

publish_mediation_report;


%% MULTIVARIATE MEDIATION
%--------------------------------------------------------------------------
% based on
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
% and personal communication with Martin Lindquist and Xiaochun Han

X_c = num2cell(X); % convert double to cell array

for i = 1:size(contrastnames,2)
    
    pdm_mediationdir_contrast = fullfile(ctq_pdm_mediationdir, contrastnames{i});
    if ~isfolder(pdm_mediationdir_contrast)
        mkdir(pdm_mediationdir_contrast);
    end
    cd(pdm_mediationdir_contrast);
    
    Y{i} = Ydat(:,i);
    Y_c{i} = num2cell(Y{i});
    
    M_brain{i} = DAT.imgs{i}(idx_ctq);
    names = M_brain{i};
    
    mask = which('gray_matter_mask.nii');
    
    for j=1:size(X_c,1)
        dat = fmri_data(names{j},mask); 
        m{j} = dat.dat; 
    end
    
    pdm = multivariateMediation(X_c,Y_c{i},m,'B',20,'svd','plots');
    pdm = multivariateMediation(pdm,'nPDM',2);
    pdm = multivariateMediation(pdm,'noPDMestimation','bootPDM',1:2,'bootjPDM','Bsamp',100);%'save2file','PDMresults.mat');
    
    dat_mask = fmri_data(mask,mask,'noverbose');
    [obj,figh] = plotPDM(pdm,dat_mask);
    
    % visualize PDMs (manually)
    % threshold PDM1 at fdr threshold from
    % pdm.pThreshold(1))
    dat.dat = pdm.Wfull{1}.*(pdm.boot.p{1}<pdm.pThreshold(1));

    % display thresholded PDM1
    figure;
    canlab_results_fmridisplay(dat,'montagetype','full');
    orthviews(dat);
    r1 = region(dat);
    table(dat);
    [r1, region_table1, table_legend_text1] = autolabel_regions_using_atlas(r1);
    
    % and PDM2
    % threshold PDM2 at fdr threshold from
    % pdm.pThreshold(2))
    dat2=dat;
    dat2.dat = pdm.Wfull{2}.*(pdm.boot.p{2}<pdm.pThreshold(2));

    % display thresholded PDM
    figure;
    canlab_results_fmridisplay(dat2,'montagetype','full');
    orthviews(dat2);
    r2 = region(dat2);
    table(dat2);
    [r2, region_table2, table_legend_text2] = autolabel_regions_using_atlas(r2);
    
    cd(ctq_pdm_mediationdir);
    
end