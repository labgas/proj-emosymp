% This script runs single-level mediation models for proj-emosymp 
% (Van den Houte, Bogaerts, et al, in prep), with the following variables
% X: group (patients versus controls)
% M: NPS response (including subregions) for each contrast
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
% This script is adapted from the following CANlab script:
%
% C:\Users\lukas\Dropbox
% (Personal)\V_visceral_common\projects\anticipation_visceral\scripts\
% s8_mediation_model_FAPS.m by @pkragel
% see also mediation walkthrough on canlab.github.io
%
% Contact @lukasvo76 if you do not have access to the original script
% and you need it
%
% @lukasvo76 @Dartmouth, March 2021
% 
% WORK IN PROGRESS
% FIGURES, TITLES, AND STRUCTURE OF HTML OUTPUT CAN BE IMPROVED
% LOOPING OVER SUBPATTERNS/REGIONS CAN BE IMPLEMENTED


%% DEFINE PATHS AND VARIABLES TO BE USED IN ALL ANALYSES
% -------------------------------------------------------------------------
addpath(genpath('C:\Users\lukas\Documents\GitHub\MediationToolbox')); % add CANlab mediation toolbox to your path
idx = ~isnan(DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu); % index for subjects with missing behavioral data, which we want to exclude
X = DAT.BETWEENPERSON.group(idx); % define predictor
Ydat = [DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_neu(idx), DAT.BEHAVIOR.behavioral_data_table.symptoms_neg_pos(idx), DAT.BEHAVIOR.behavioral_data_table.symptoms_pos_neu(idx)]; % define matrix of outcomes (subjects * contrasts)


%% NPS
% -------------------------------------------------------------------------
contrastnames = DAT.SIG_contrasts.raw.dotproduct.conditionnames; % get names of contrasts

for i = 1:size(contrastnames,2) % loop over contrasts
    
    % define outcome for contrast i
    Y{i} = Ydat(:,i);
    
    % define mediators
    M{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPS(idx,i));
    Mpos{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSpos(idx,i));
    Mneg{i} = table2array(DAT.SIG_contrasts.raw.dotproduct.NPSneg(idx,i));
    
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
contrastnames = DAT.SIG_contrasts.raw.dotproduct.conditionnames;

for i = 1:size(contrastnames,2) 
    
    Y{i} = Ydat(:,i);
    
    MrIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,2);
    MrdpIns{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,6);
    MrS2_Op{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,7);
    MdACC{i} = DAT.NPSsubregions.npspos_by_region_contrasts{i}(idx,8);
    
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
    
end