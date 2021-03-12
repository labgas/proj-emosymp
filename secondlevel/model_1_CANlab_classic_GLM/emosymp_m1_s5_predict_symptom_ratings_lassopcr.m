% This script runs a 5-fold cross-validated lasso-PCR model to predict
% somatic symptom ratings from brain responses to positive, negative, and
% neutral IAPS pictures (3 condition images from classic GLM per subject)
% for proj-emosymp (Van den Houte, Bogaerts, et al, in prep).
%
% It requires all study-specific and core prep_ scripts from proj-emosymp 
% (https://github.com/labgas/proj-emosymp) and CANlab_help examples (LaBGAS
% fork, https://github.com/labgas/CANlab_help_examples) 
% to be run first, or load your saved .mat files if you ran these
% scripts earlier!
%
% Options for masking and scaling are set in
% a2_emosymp_m1_s2_set_default_options
%
% This script is adapted from the following CANlab scripts
%
% 1. C:\Users\lukas\Dropbox
%   (Personal)\V_visceral_common\projects\anticipation_visceral\scripts\
%   s1_predict_anxiety_ratings_lassopcr.m by @pkragel
% 2. C:\Users\lukas\Dropbox (Dartmouth College)\CANlab walkthroughs\
%   CANlab_single_trials_demo_14\CANlab_single_trials_demo_14.m by
%   @bogpetre and @lukasvo76 - see also walkthrough on canlab.github.io
% 3. C:\Users\lukas\Dropbox (Dartmouth College)\CANlab walkthroughs\
%   CANlab_MVPA_between_within\CANlab_MVPA_between_within.m by @bogpetre
%   and @lukasvo76 - see also walkthrough on canlab.github.io
%
% Contact @lukasvo76 if you do not have access to these original scripts
% (and their helper functions) and you need them
%
% @lukasvo76 @Dartmouth, March 2021
% 
% WORK IN PROGRESS
% FIGURES, TITLES, AND STRUCTURE OF HTML OUTPUT NEEDS TO BE IMPROVED!


%% CHECK OPTIONS FROM A_2_SET_DEFAULT_OPTIONS
%--------------------------------------------------------------------------
options_needed = {'dosavepcrstats', 'dobootstrap', 'boot_n', 'myscaling_pcr'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {true false 1000 'raw'};          % defaults if we cannot find info in a2_set_default_options at all 

plugin_get_options_for_analysis_script


%% CHECK SPIDER TOOLBOX
% -------------------------------------------------------------------------

spath = which('use_spider.m');
if isempty(spath)
    disp('Warning: spider toolbox not found on path; prediction may break')
end


%% GET MASK
% -------------------------------------------------------------------------

if exist('maskname_pcr', 'var') && ~isempty(maskname_pcr)
    pcrmask = fmri_data(maskname_pcr, 'noverbose');
end


%% SET UP AND TRAIN PREDICTION MODEL ON CONDITIONS
% -------------------------------------------------------------------------

c = size(DAT.conditions, 2);
wh = [1:c];

behdat = DAT.BEHAVIOR.behavioral_data_table;
norating_idx = ~isnan(behdat.symptoms_neg);
behdat = behdat(norating_idx,:);
group_idx = DAT.BETWEENPERSON.group(norating_idx);
subj_idx = behdat.id;

if dobootstrap, pcrtime = tic; end

% concatenate data objects for each condition into one
switch myscaling_pcr
    case 'raw'
        scaling_string = 'no_scaling';
        for i = 1:c
            DATA_OBJ_rating{i} = DATA_OBJ{i}.get_wh_image(norating_idx');
        end
        [cat_obj, condition_idx] = cat(DATA_OBJ_rating{wh});
    case 'scaled'
        scaling_string = 'scaling_z_score_conditions';
        for i = 1:c
            DATA_OBJsc_rating{i} = DATA_OBJsc{i}.get_wh_image(norating_idx');
        end
        [cat_obj, condition_idx] = cat(DATA_OBJsc_rating{wh});
    otherwise
        error('myscaling must be ''raw'' or ''scaled''');
end

cat_obj = remove_empty(cat_obj);

% apply mask if specified in a2_ script
if exist('pcrmask', 'var')        
        disp('Masking data')
        cat_obj = apply_mask(cat_obj, pcrmask);        
else
        disp('No mask found; using full existing image data');     
end

% add outcome variable to fmri_data_st object
cat_obj.Y = [behdat.symptoms_neg; behdat.symptoms_neu; behdat.symptoms_pos];
cat_obj.Y = log(cat_obj.Y); % taking the natural log here because of right tail

% define holdout sets, with the following rules
% 1. balanced over groups
% 2. stratified for conditions, i.e. leave entire subject out
% this combination can easily be achieved with @bogpetre's cvpartition2
% object - check help for cvpartition and cvpartition2 for more info
nfolds = 5; % @lukasvo76: hardcoded for now, but can easily be adapted
groupc = repmat(group_idx,c);
groupc = groupc(:,1);
subjc = repmat(subj_idx,c);
subjc = subjc(:,1);
cvpart = cvpartition2(groupc,'KFOLD',nfolds,'Stratify',subjc);
[I,J] = find([cvpart.test(1),cvpart.test(2), cvpart.test(3), cvpart.test(4), cvpart.test(5)]);
holdout_set = sortrows([I,J]);
holdout_set = holdout_set(:,2);

% in case you only want to include certain conditions, not all of them
% cat_obj=cat_obj.get_wh_image(condition_idx' < 3);
% holdout_set=holdout_set(condition_idx < 3);

% train model
if dobootstrap
    [pcr_cverr, pcr_stats] = predict(cat_obj, 'algorithm_name', 'cv_lassopcr', 'nfolds', holdout_set, 'bootsamples', boot_n);
else
    [pcr_cverr, pcr_stats] = predict(cat_obj, 'algorithm_name', 'cv_lassopcr', 'nfolds', holdout_set);
end

% print and plot results
fprintf('PCR r = %0.3f\n', corr(pcr_stats.yfit, cat_obj.Y));

figure
line_plot_multisubject(cat_obj.Y, pcr_stats.yfit, 'subjid', condition_idx');
xlabel({'Observed Symptoms','(stim level average)'}); ylabel({'PCR Estimated Symptoms','(cross validated)'});
drawnow; snapnow;
% @lukasvo76: note that within- and between-person r in output reflects
% within- and between-conditions, since we put condition_idx as subject
% identifier to plot line per condition

figure
pcr_stats.weight_obj.montage;
drawnow;snapnow;

% save stats
if dosavepcrstats
    pcr_stats_results_conditions = pcr_stats;
    
        if exist('pcrmask', 'var')

            pcr_stats_results_conditions.mask = pcrmask;
            pcr_stats_results_conditions.maskname = maskname_pcr;

        end
        
    savefilenamedata = fullfile(resultsdir,['pcr_stats_results_conditions_',scaling_string,'.mat']);

    save(savefilenamedata, 'pcr_stats_results_conditions', '-v7.3');
    printhdr('Saved pcr_stats_results for conditions');
    
end


if dobootstrap, disp('Cumulative run time:'), toc(svmtime);end


%% SET UP AND TRAIN PREDICTION MODEL ON CONTRASTS
% -------------------------------------------------------------------------
clear cat_obj

kc = size(DAT.contrastnames, 2);
kwh = [1:kc];

pcr_stats_results_contrasts = cell(1, kc);

behdat = DAT.BEHAVIOR.behavioral_data_table;
norating_idx = ~isnan(behdat.symptoms_neg_neu);
behdat = behdat(norating_idx,:);
group_idx = DAT.BETWEENPERSON.group(norating_idx);
subj_idx = behdat.id;
symptom_ratings = [behdat.symptoms_neg_neu, behdat.symptoms_neg_pos, behdat.symptoms_pos_neu];

for i = 1:kc

if dobootstrap, pcrtime = tic; end
    
switch myscaling_pcr
    case 'raw'
        cat_obj{i} = DATA_OBJ_CON{i}.get_wh_image(norating_idx');
    case 'scaled'
        cat_obj{i} = DATA_OBJ_CONsc{i}.get_wh_image(norating_idx');
    case 'scaled_contrasts'
        cat_obj{i} = DATA_OBJ_CONscc{i}.get_wh_image(norating_idx');
    otherwise
        error('myscaling must be ''raw'' or ''scaled'' or ''scaled_contrasts''');
    end

% apply mask if specified in a2_ script
if exist('pcrmask', 'var')        
        disp('Masking data')
        cat_obj{i} = apply_mask(cat_obj{i}, pcrmask);        
else
        disp('No mask found; using full existing image data');     
end

% add outcome variable to fmri_data_st object
cat_obj{i}.Y = symptom_ratings(:,i);
cat_obj{i}.Y = log(1 + (cat_obj{i}.Y - min(cat_obj{i}.Y)));

% define holdout sets, balanced over groups
% this combination can easily be achieved with @bogpetre's cvpartition2
% object - check help for cvpartition and cvpartition2 for more info
nfolds = 5; % @lukasvo76: hardcoded for now, but can easily be adapted
cvpart = cvpartition(group_idx,'KFOLD',nfolds);
[I,J] = find([cvpart.test(1),cvpart.test(2), cvpart.test(3), cvpart.test(4), cvpart.test(5)]);
holdout_set = sortrows([I,J]);
holdout_set = holdout_set(:,2);

% code in case you only want to include certain conditions, not all of them
% cat_obj=cat_obj.get_wh_image(condition_idx' < 3);
% holdout_set=holdout_set(condition_idx < 3);

% train model
if dobootstrap
    [pcr_cverr, pcr_stats] = predict(cat_obj{i}, 'algorithm_name', 'cv_lassopcr', 'nfolds', holdout_set, 'bootsamples', boot_n);
else
    [pcr_cverr, pcr_stats] = predict(cat_obj{i}, 'algorithm_name', 'cv_lassopcr', 'nfolds', holdout_set);
end

% print and plot results
fprintf('PCR r = %0.3f\n', corr(pcr_stats.yfit, cat_obj{i}.Y));

figure
scatter(cat_obj{i}.Y, pcr_stats.yfit);
xlabel({'Observed Symptoms','(stim level average)'}); ylabel({'PCR Estimated Symptoms','(cross validated)'});
drawnow;snapnow;

figure
pcr_stats.weight_obj.montage;
drawnow;snapnow;

% store stats in cell array
pcr_stats_results_contrasts{i} = pcr_stats;
    
        if exist('pcrmask', 'var')

            pcr_stats_results_contrasts{i}.mask = pcrmask;
            pcr_stats_results_contrasts{i}.maskname = maskname_pcr;

        end
        
if dobootstrap, disp('Cumulative run time:'), toc(svmtime);end

end

% save stats
if dosavepcrstats
        
    savefilenamedata = fullfile(resultsdir,['pcr_stats_results_contrasts_',scaling_string,'.mat']);

    save(savefilenamedata, 'pcr_stats_results_contrasts', '-v7.3');
    printhdr('Saved pcr_stats_results for contrasts');
    
end