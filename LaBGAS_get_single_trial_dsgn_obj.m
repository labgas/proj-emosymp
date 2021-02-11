%% LaBGAS_get_firstlvl_dsgn_obj
%
% This script contains a function that defines a number of fields of 
% the CANlab style first level DSGN structure array
% for Maaike's emosymp study on IAPS picture-induced bodily symptoms.
% This function is used in LaBGAS_1_spm_fit_firstlvl_models.m to
% run first level analysis using CANlab tools
%
% THIS SCRIPT IS IDENTICAL TO LaBGAS_get_firstlvl_dsgn_obj.m except for
% 1) single trial settings in DSGN.singletrials!
% 2) names of directories in DSGN.modeldir and DSGN.modelingfilesdir 
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
% date:   October, 2020
%
%__________________________________________________________________________
% @(#)% LaBGAS_get_firstlvl_dsgn_obj.m         v1.1        
% last modified: 2021/02/11
%
%% function code
function DSGN = LaBGAS_get_single_trial_dsgn_obj()
    % INPUT
    % required fields
    DSGN.metadata = "emosymp study classic GLM first level analysis using CANlab tools"; % field for annotation with study info, or whatever you like
    DSGN.modeldir = 'C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\firstlevel\model_2_CANlab_single_trial_glm'; % directory where you want to write first level results
    DSGN.subjects = {}; % sets up empty cell array field for subjects in structural array DSGN
    fnames = dir('C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\derivatives\fmriprep\sub*'); % get subject names from root directory with subject data
    idx=[fnames.isdir]'; % duplicate names of subjects in fmriprep dir because of .html files with same subject names as folders
    fnames=fnames(idx);
    for i = 1:size(fnames,1)
        this_f = fnames(i);
            runs = dir([this_f.folder, '\', this_f.name, '\func\sub*.tsv']);
            if length(runs) > 1 % only include subject if it has at least two runs - this should be the case in all subjects for this study - check trouble_in_paradise.xlsx
                DSGN.subjects = [DSGN.subjects, [this_f.folder, '\', this_f.name]];
            end
    end
    DSGN.funcnames = {'func\run-1\s6*.nii',...
        'func\run-2\s6*.nii',...
        'func\run-3\s6*.nii',...
        'func\run-4\s6*.nii',...
        'func\run-5\s6*.nii',...
        'func\run-6\s6*.nii'}; % cell array of subject directories (absolute paths)	
    % optional fields
    DSGN.concatenation = {}; % default: none; cell array of arrays of runs to concatenate; see documentation for when to concatenate, and how it works exactly
    DSGN.allowmissingfunc = true; % default; true will prevent erroring out when functional file is missing
    DSGN.customrunintercepts = {}; % default: none; will only work if DSGN.concatenation is specified; cell array of vectors specifying custom intercepts
    
    % PARAMETERS
    DSGN.tr = 3; % repetition time (TR) in seconds
    DSGN.hpf = 180; % high pass filter in seconds; SPM default is 128, CANlab default is 180 since the brain response to pain stimuli last long and variance may be lost at shorter lengths, use scn_spm_design_check output for diagnostics; STUDY-SPECIFIC
    DSGN.fmri_t = 50; % microtime resolution - t=number of slices since we did slice timing; spm (and CANlab) default 16 can be kept for multiband w/o slice timing
    DSGN.fmri_t0 = 25; % microtime onset - reference slice used in slice timing correction; spm (and CANlab) default 1 can be kept for multiband w/o slice timing
    
    % MODELING
    % required fields: cell array (one cell per session) of cell arrays (one cell per condition) of MAT-file names; if only one session is specified, it will be applied to all sessions
    c=0;
    c=c+1;DSGN.conditions{1}{c}='negative';
    c=c+1;DSGN.conditions{1}{c}='neutral';
    c=c+1;DSGN.conditions{1}{c}='positive';
    c=c+1;DSGN.conditions{1}{c}='scoring';
    % optional fields
%     DSGN.pmods = {{}}; % cell array (one cell per session) of cell arrays (one cell per condition) of cell arrays (one cell per modulator) of MAT-file names
%     DSGN.convolution; default hrf.derivs = [0 0]; structure specifying the convolution to use for conditions different fields required depending on convolution type; 
%     DSGN.ar1 = false; % autoregressive AR(1) to model serial correlations; SPM default is true, CANlab default is false, Tor recommends turning autocorrelation off, because this algorithm pools across the whole brain, and does not perform well in some situations; if you are performing a group analysis, the autocorrelation problem is not as concerning
    DSGN.notimemod = true; % default: false; if true, turn off time modulation of conditions, i.e. when you do not expect linear trends over time
    DSGN.singletrials{1} = {1 1 1 0}; % a cell array (1 cell per session) of cell arrays (1 cell per condition) of (corresponding to DSGN.conditions) of true/false values indicating whether to convert specified condition to set of single trial conditions
%     DSGN.singletrialsall = false; % default: false; if true, set DSGN.singletrials to true for all conditions
    DSGN.modelingfilesdir = 'model_2_CANlab_single_trial_glm'; % name of subfolder which will be created within directory containing functional files where .mat files containing fields of DSGN structure will be saved
%     DSGN.allowemptycond = false; % default:false; if true, allow empty conditions
%     DSGN.allowmissingcondfiles = false; % default:false; if true, throw warning instead of error when no file(s) are found corresponding to a MAT-file name/wildcard
    DSGN.multireg = 'noise_regs'; % specify name for matfile with noise parameters you want to save
    
    % CONTRASTS
    % required fields
    % cell array (one cell per contrast) of contrast definitions
    c=0;
    c=c+1;DSGN.contrasts{c} = {{'negative'}};
    c=c+1;DSGN.contrasts{c} = {{'neutral'}};
    c=c+1;DSGN.contrasts{c} = {{'positive'}};
    c=c+1;DSGN.contrasts{c} = {{'negative'} {'neutral'}};
    c=c+1;DSGN.contrasts{c} = {{'negative'} {'positive'}};
    c=c+1;DSGN.contrasts{c} = {{'positive'} {'neutral'}};
    % optional fields - not needed for standard contrasts
%     DSGN.contrastnames{1} = 'negative'; 
%     DSGN.contrastweights{1} = [1]; 
%     DSGN.contrastnames{2} = 'neutral'; 
%     DSGN.contrastweights{2} = [0 1]; 
%     DSGN.contrastnames{3} = 'positive'; 
%     DSGN.contrastweights{3} = [0 0 1];
%     DSGN.contrastnames{4} = 'negative versus neutral'; 
%     DSGN.contrastweights{4} = [1 -1];
%     DSGN.contrastnames{5} = 'negative versus positive'; 
%     DSGN.contrastweights{5} = [1 0 -1];
%     DSGN.contrastnames{6} = 'positive versus neutral'; 
%     DSGN.contrastweights{6} = [0 -1 1];
end