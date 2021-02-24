% Always run this first before you run other scripts.
%
% NOTES:
% - standard folders and variable names are created by these scripts
%
% - in "prep_" scripts: 
%   image names, conditions, contrasts, colors, global gray/white/CSF
%   values are saved automatically in a DAT structure
% 
% - extracted fmri_data objects are saved in DATA_OBJ variables
% - contrasts are estimated and saved in DATA_OBJ_CON variables
%
% - files with these variables are saved and loaded automatically when you
%   run the scripts
%   meta-data saved in image_names_and_setup.mat
%   image data saved in data_objects.mat
%
% - you only need to run the prep_ scripts once.  After that, use 
%   b_reload_saved_matfiles.m to re-load saved files
% 
% - when all scripts working properly, run z_batch_publish_analyses.m
%   to create html report.  customize by editing z_batch_list_to_publish.m
%
% - saved in results folder:
%   figures
%   html report with figures and stats, in "published_output"
%
% LUKASVO76's NOTES
% - this script is copied and edited from 
%   https://github.com/canlab/CANlab_help_examples/tree/master/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify
% - as described in the documentation in https://canlab.github.io/batch/, 
%   you can first type a_set_up_new_analysis_folder_and_scripts 
%   in your Matlab command window and run it to automatically 
%   set up folder structure according to CANlab conventions and
%   copy in some of the CANlab template scripts, including the present one,
%   automatically; an alternative is to copy and edit that script
%   to adapt your folder structure for example according to BIDS convention
% - to check whether the correct dependencies are on your Matlab path, you
%   may want to run a2_second_level_toolbox_check_dependencies from your
%   Matlab command line as is from the CANlab_help_examples repo, no need
%   for study-specific edits!
% - if you wish to change the default options for some of the prep scripts,
%   you should copy the script a2_set_default_options which is called by
%   this script to your local scripts directory, rename and adapt what you want
%   DO NOT CHANGE THE VERSION OF THAT SCRIPT IN THE CANLAB OR LABGAS GITHUB FOLDER 
%   WHERE IT LIVES ORIGINALLY, at least not without branching or forking!

% Set base directory
% --------------------------------------------------------

% lukasvo76: Base directory for your second level model

basedir = 'C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\secondlevel\model_1_CANlab_classic_GLM';

% Set user options
% --------------------------------------------------------

a2_emosymp_m1_s2_set_default_options % lukasvo76: referring to renamed study- and model-specific script here, as I want to customize the default option - see my note in documentation above

% Set up paths
% --------------------------------------------------------

cd(basedir)

datadir = fullfile(basedir, 'data');
resultsdir = fullfile(basedir, 'results');
scriptsdir = 'C:\Users\lukas\Documents\GitHub\proj-emosymp\secondlevel\model_1_CANlab_classic_GLM'; % lukasvo76: contrary to CANlab structure, I want my second level scripts to live in my Github local working dir for this study
figsavedir = fullfile(resultsdir, 'figures');
notesdir = fullfile(resultsdir, 'notes');
canlabhelpexamplesdir = 'C:\Users\lukas\Documents\GitHub\CANlab_help_examples_LaBGAS_fork\CANlab_help_examples'; % lukasvo76: we want the LaBGAS fork of CANlab_help_examples here!

addpath(scriptsdir)
addpath(genpath(canlabhelpexamplesdir))

if ~exist(resultsdir, 'dir'), mkdir(resultsdir); end
if ~exist(figsavedir, 'dir'), mkdir(figsavedir); end

% Display helper functions: Called by later scripts
% --------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);
