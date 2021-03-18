%% LaBGAS_extract_pmods_phenotype_trial_tsv.m
%
% This script will extract symptom ratings as parametric modulators 
% from the phenotype_trial.tsv file for each condition of the emosymp study 
% 
% DEPENDENCIES
% BIDS data organization
%
% INPUTS 
% phenotype_trial.tsv file
%
% OUTPUT
% noise_regs & onsets files that can be loaded into CANlab DSGN structure
% (RECOMMENDED)
% or directly into SPM first level batch (for the latter, use
% LaBGAS_first_level_batch_fMRIprep_conf)
%
% IMPORTANT NOTE!
% for future studies, we may want to integrate this script with:
% 1. LaBGAS_BIDS_emosymp_fMRI - last section writing the events.tsv files
% 2. LaBGAS_extract_confound_reg_fMRIprep
% to create the following files in one script, 
% and put them in the right place in the CANlab directory structure:
% 1. events.tsv (including ratings/pmods) from logfiles
% 2. phenotype_trial.tsv (for behavioral analysis purposes in SAS/Matlab)
%    from logfiles
% 3. noise_reg, onset/duration, and pmod files (as needed for CANlab first
%    level batch tools) from the above .tsv files
%
% let @lukasvo76 know if you want to implement this for your study, 
% so we can work on it!
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   March, 2021
%
%__________________________________________________________________________
% @(#)% LaBGAS_extract_pmods_phenotype_trial_tsv.m         v1.0        
% last modified: 2021/03/18


%% define dirs
rootdir='C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS';
derivdir=fullfile(rootdir,'derivatives\fmriprep');
phenodir=fullfile(rootdir,'phenotype');


%% load phenotype_trial.tsv into Matlab as table, define subjects, and perform sanity check
cd(phenodir);
pheno_trial = readtable('Phenotype_trial.tsv','FileType', 'text', 'Delimiter', 'tab');
idx_subjincluded = logical(pheno_trial.participant_included);
idx_runincluded = logical(pheno_trial.run_included);
idx_norating = ~isnan(pheno_trial.state_symptoms);
idx = idx_subjincluded & idx_runincluded & idx_norating;
pheno_trial = pheno_trial(idx,:);
pheno_trial.subject_id = categorical(pheno_trial.subject_id);
pheno_trial.diagnosis = categorical(pheno_trial.diagnosis);
pheno_trial.trial_valence = categorical(pheno_trial.trial_valence);
subjs=unique(pheno_trial.subject_id);
% sanity check on subject definition
pheno_sub = readtable('Phenotype_subject.tsv','FileType','text','Delimiter','tab');
idx_subjincluded_sub = logical(pheno_sub.included);
pheno_sub = pheno_sub(idx_subjincluded_sub,:);
idx_norating_sub = ~isnan(pheno_sub.symptoms_neg) & ~isnan(pheno_sub.symptoms_neu) & ~isnan(pheno_sub.symptoms_pos);
pheno_sub = pheno_sub(idx_norating_sub,:);
cd(derivdir);
subs=dir('sub-*');
id=[subs.isdir]';
subs={subs(id).name}';
subs=subs(idx_norating_sub);
if subjs == subs
    disp('subject names in phenotype_trial.tsv match names of subject dirs in derivdir - proceeding');
else
    error('subject names in phenotype_trial.tsv do not match names of subject dirs in derivdir - check before proceeding')
end
subjs=cellstr(subjs);


%% loop over subjects
for sub=1:size(subjs,1)
    % define subjectdir & rundirs
    subjderivdir=fullfile(derivdir,subjs{sub},'func');
    cd(subjderivdir);
    subj_rundirs=ls('run-*');
    
    % get subject-level data
    subj_name = subjs{sub};
    subj_pheno = pheno_trial(pheno_trial.subject_id==subj_name,:);
    subj_pheno.state_symptoms = subj_pheno.state_symptoms - mean(subj_pheno.state_symptoms); % mean center your pmod at the subject level; comment out if you don't want to do that
    subj_runs = unique(subj_pheno.run)';
    
    % loop over runs
    for run=subj_runs
        cd(subjderivdir);
        run_pheno = subj_pheno(subj_pheno.run==run,:);
        cd(subj_rundirs(run,:));
        % save pmods file as .mat file
        filename = strcat('pmods_',subj_name,'_',subj_rundirs(run,:));
        save(filename,'run_pheno');
    end
end