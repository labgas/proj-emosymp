%% LaBGAS_get_firstlvl_dsgn_obj_pmods
%
% This script contains a function that defines a number of fields of 
% the CANlab style first level DSGN structure array
% for Maaike's emosymp study on IAPS picture-induced bodily symptoms, 
% including parametric modulators.
% This function is used in LaBGAS_1_spm_fit_firstlvl_models.m to
% run first level analysis using CANlab tools
%
% IMPORTANT NOTE: function and script are study-specific!
%
% See canlab_glm_single_subject('dsgninfo')
% OR Github\CanlabCore\CanlabCore\GLM_Batch_tools\canlab_glm_dsgninfo.txt 
% for details on first level analysis using CANlab tools, 
% including options and defaults
%
% This script is adapted by @lukasvo76 from the scripts
% 1) get_firslvl_dsgn_obj.m and get_single_trial_dsgn_obj.m by @bogpetre on
% Google Drive\CANlab\CANLAB Lab Member Documents\GLM_batch_tools\
% bogdan_paingen\classic_glm_contrasts
% 2) MPA2_set_design_model1_blanca.m by @martaceko on
% Google Drive\CANlab\CANLAB Lab Member Documents\GLM_batch_tools\
% Marta_MPA2\MPA2code_1stlevel\Code
% contact @lukasvo76 if you need those original scripts
% 
% DEPENDENCIES ON YOUR MATLAB PATH
% a) SPM12
% b) CANlab tools cloned from Github (see canlab.github.io)
% 
% INPUTS 
% none - you need to adapt the function script below to your study
%
% OUTPUT
% CANlab style first level DSGN structure array in your Matlab workspace, 
% to be used by LaBGAS_1_spm_fit_firstlvl_models.m
%
%__________________________________________________________________________
%
% authors: 
% lukas.van.oudenhove@dartmouth.edu, lukas.vanoudenhove@kuleuven.be
% bogdan.petre@dartmouth.edu,
% marta.ceko@colorado.edu
%
% date:   March, 2021
%
%__________________________________________________________________________
% @(#)% LaBGAS_get_firstlvl_dsgn_obj_pmods.m         v1.0        
% last modified: 2021/03/16
%
%% function code
function DSGN = LaBGAS_get_firstlvl_dsgn_obj_pmods()
    % INPUT
    % required fields
    DSGN.metadata = "emosymp study classic GLM first level analysis using CANlab tools"; % field for annotation with study info, or whatever you like
    DSGN.modeldir = 'C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\firstlevel\model_3_CANlab_classic_glm_pmods'; % directory where you want to write first level results
    DSGN.subjects = {}; % sets up empty cell array field for subjects in structural array DSGN
    fnames = dir('C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\derivatives\fmriprep\sub*'); % get subject names from root directory with subject data
    idx=[fnames.isdir]'; % duplicate names of subjects in fmriprep dir because of .html files with same subject names as folders
    fnames=fnames(idx);
    for i = 1:size(fnames,1)
        this_f = fnames(i);
            runs = dir([this_f.folder, '\', this_f.name, '\func\sub*.tsv']);
            if length(runs) > 1 % only include subject if it has at least two runs - this should be the case in all subjects for this study - check trouble_in_paradise.xlsx
                DSGN.subjects = [DSGN.subjects, [this_f.folder, '\', this_f.name]]; % cell array of subject directories (absolute paths)
            end
    end
    pheno_sub = readtable('C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\phenotype\Phenotype_subject.tsv','FileType','text','Delimiter','tab');
    idx_subjincluded_sub = logical(pheno_sub.included);
    pheno_sub = pheno_sub(idx_subjincluded_sub,:);
    idx_norating_sub = ~isnan(pheno_sub.symptoms_neg) & ~isnan(pheno_sub.symptoms_neu) & ~isnan(pheno_sub.symptoms_pos);
    DSGN.subjects = DSGN.subjects(idx_norating_sub'); % get rid of subjects who have missing pmods
    DSGN.funcnames = {'func\run-1\s6*.nii',...
        'func\run-2\s6*.nii',...
        'func\run-3\s6*.nii',...
        'func\run-4\s6*.nii',...
        'func\run-5\s6*.nii',...
        'func\run-6\s6*.nii'}; % cell array (one cell per session) of paths to functional files, relative to absolute path specific in DSGN.subjects
    % optional fields
    DSGN.concatenation = {[1:6]}; % default: none; cell array of arrays of runs to concatenate; see documentation for when to concatenate, and how it works exactly; in this case we concatenate for several of the reasons outlined in https://docs.google.com/document/d/1pfQlFnumT1ZjEB3qSHWA9lZ05Zd6pDqxwUy3jSNOhnY/edit#heading=h.i5hmh8orkhrs
    DSGN.allowmissingfunc = true; % default; true will prevent erroring out when functional file is missing
    DSGN.customrunintercepts = {}; % default: none; will only work if DSGN.concatenation is specified; cell array of vectors specifying custom intercepts; not needed in this case, concatenation automatically adds run intercept for each run!
    
    % PARAMETERS
    DSGN.tr = 3; % repetition time (TR) in seconds
    DSGN.hpf = 180; % high pass filter in seconds; SPM default is 128, CANlab default is 180 since the brain response to pain stimuli last long and variance may be lost at shorter lengths, use scn_spm_design_check output for diagnostics; STUDY-SPECIFIC
    DSGN.fmri_t = 50; % microtime resolution - t=number of slices since we did slice timing; spm (and CANlab) default 16 can be kept for multiband w/o slice timing
    DSGN.fmri_t0 = 25; % microtime onset - reference slice used in slice timing correction; spm (and CANlab) default 1 can be kept for multiband w/o slice timing
    
    % MODELING
    % required fields: cell array (one cell per session) of cell arrays (one cell per condition) of MAT-file names; if only one session is specified, it will be applied to all sessions
    c=0;
    c=c+1;DSGN.conditions{c}={'negative' 'neutral' 'positive' 'scoring'};
    c=c+1;DSGN.conditions{c}={'negative' 'neutral' 'positive' 'scoring'};
    c=c+1;DSGN.conditions{c}={'negative' 'neutral' 'positive' 'scoring'};
    c=c+1;DSGN.conditions{c}={'negative' 'neutral' 'positive' 'scoring'};
    c=c+1;DSGN.conditions{c}={'negative' 'neutral' 'positive' 'scoring'};
    c=c+1;DSGN.conditions{c}={'negative' 'neutral' 'positive' 'scoring'};
    % optional fields
    c=0;
    c=c+1;DSGN.pmods{c}={'symptoms_neg' 'symptoms_neu' 'symptoms_pos'}; % cell array (one cell per session) of cell arrays (one cell per condition) of cell arrays (one cell per modulator) of MAT-file names - corresponds to DSGN.conditions, except for the conditions you don't want to add a parametric modulator for
    c=c+1;DSGN.pmods{c}={'symptoms_neg' 'symptoms_neu' 'symptoms_pos'};
    c=c+1;DSGN.pmods{c}={'symptoms_neg' 'symptoms_neu' 'symptoms_pos'};
    c=c+1;DSGN.pmods{c}={'symptoms_neg' 'symptoms_neu' 'symptoms_pos'};
    c=c+1;DSGN.pmods{c}={'symptoms_neg' 'symptoms_neu' 'symptoms_pos'};
    c=c+1;DSGN.pmods{c}={'symptoms_neg' 'symptoms_neu' 'symptoms_pos'};
%     DSGN.convolution; default hrf.derivs = [0 0]; structure specifying the convolution to use for conditions different fields required depending on convolution type; 
%     DSGN.ar1 = false; % autoregressive AR(1) to model serial correlations; SPM default is true, CANlab default is false, Tor recommends turning autocorrelation off, because this algorithm pools across the whole brain, and does not perform well in some situations; if you are performing a group analysis, the autocorrelation problem is not as concerning
    DSGN.notimemod = true; % default: false; if true, turn off time modulation of conditions, i.e. when you do not expect linear trends over time
%     DSGN.singletrials = {{}}; % a cell array (1 cell per session) of cell arrays (1 cell per condition) of (corresponding to DSGN.conditions) of true/false values indicating whether to convert specified condition to set of single trial conditions
%     DSGN.singletrialsall = false; % default: false; if true, set DSGN.singletrials to true for all conditions
    DSGN.modelingfilesdir = 'model_3_CANlab_classic_glm_pmods'; % name of subfolder which will be created within directory containing functional files where .mat files containing fields of DSGN structure will be saved
%     DSGN.allowemptycond = false; % default:false; if true, allow empty conditions
%     DSGN.allowmissingcondfiles = false; % default:false; if true, throw warning instead of error when no file(s) are found corresponding to a MAT-file name/wildcard
    DSGN.multireg = 'noise_regs'; % specify name for matfile with noise parameters you want to save
    
    % CONTRASTS
    % optional fields for flexible definition of contrasts - for more info
    % help canlab_spm_contrast_job_luka
    % https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_example_DSGN_setup.txt
    DSGN.regmatching = 'regexp'; % regular expression mode to match keywords you provide in cell arrays below with beta regressor names stored in the SPM.Vbeta.descrip field of your first level SPM.mat file
%     DSGN.defaultsuffix = '\*bf\(1\)$'; % adds this suffix to each keyword your provide by default
    % required fields
    % cell array (one cell per contrast) of contrast definitions
    c=0;
    c=c+1;DSGN.contrasts{c} = {{'.*negative{1}\s[^x]'}}; % this will select any beta regressor starting with "negative", followed by exactly one white space, but not followed by x - which is only the unmodulated regressors for the negative condition!
    c=c+1;DSGN.contrasts{c} = {{'.*neutral{1}\s[^x]'}};
    c=c+1;DSGN.contrasts{c} = {{'.*positive{1}\s[^x]'}};
    c=c+1;DSGN.contrasts{c} = {{'.*negative{1}\s[^x]'} {'.*neutral{1}\s[^x]'}}; % contrasts between unmodulated regressors
    c=c+1;DSGN.contrasts{c} = {{'.*negative{1}\s[^x]'} {'.*positive{1}\s[^x]'}};
    c=c+1;DSGN.contrasts{c} = {{'.*positive{1}\s[^x]'} {'.*neutral{1}\s[^x]'}};
    c=c+1;DSGN.contrasts{c} = {{'.*symptoms_neg'}}; % this will select any beta regressor with "symptoms_neg" anywhere in the name - which is only the modulated regressors for the negative condition!
    c=c+1;DSGN.contrasts{c} = {{'.*symptoms_neu'}};
    c=c+1;DSGN.contrasts{c} = {{'.*symptoms_pos'}};
    c=c+1;DSGN.contrasts{c} = {{'.*symptoms_neg'} {'.*symptoms_neu'}}; % contrasts between modulated regressors
    c=c+1;DSGN.contrasts{c} = {{'.*symptoms_neg'} {'.*symptoms_pos'}};
    c=c+1;DSGN.contrasts{c} = {{'.*symptoms_pos'} {'.*symptoms_neu'}};
    % optional fields to define custom contrast names and weights
    c=0;
    c=c+1;DSGN.contrastnames{c} = 'negative unmodulated';  
    c=c+1;DSGN.contrastnames{c} = 'neutral unmodulated';  
    c=c+1;DSGN.contrastnames{c} = 'positive unmodulated';
    c=c+1;DSGN.contrastnames{c} = 'negative unmodulated versus neutral unmodulated';  
    c=c+1;DSGN.contrastnames{c} = 'negative unmodulated versus positive unmodulated';  
    c=c+1;DSGN.contrastnames{c} = 'positive unmodulated versus neutral unmodulated';
    c=c+1;DSGN.contrastnames{c} = 'negative modulated';  
    c=c+1;DSGN.contrastnames{c} = 'neutral modulated';  
    c=c+1;DSGN.contrastnames{c} = 'positive modulated';
    c=c+1;DSGN.contrastnames{c} = 'negative modulated versus neutral modulated';  
    c=c+1;DSGN.contrastnames{c} = 'negative modulated versus positive modulated';  
    c=c+1;DSGN.contrastnames{c} = 'positive modulated versus neutral modulated';
end