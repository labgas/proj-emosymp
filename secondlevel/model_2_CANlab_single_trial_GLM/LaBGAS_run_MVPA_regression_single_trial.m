%% LaBGAS_run_MVPA_regression_single_trial
%
% This script runs MVPA regression analysis on a continuous outcome Y
% (lasso-pcr, but can easily be adapted to support vector regression or
% other methods available in CANlab's predict function) on an fmri_data_st
% object created using LaBGAS_create_single_trial_fmri_data_st_obj.m - that
% script should be run first, or the present script will load the data
% object if it is saved by the previous script.
%
% It is based on the extremely helpful tutorials on single trial analysis
% in the context of MVPA by @bogpetre @CANlab.
%
% Here are Bogdan's walkthroughs:
% https://canlab.github.io/_pages/canlab_single_trials_demo/demo_norming_comparison.html
% https://canlab.github.io/_pages/mlpcr_demo/mlpcr_demo.html (WiP)
%
% Here are two scripts Lukas adapted from these walkthroughs
% https://www.dropbox.com/sh/e17nl3ew1db1twk/AACO9QAEt6Sy3TejH-n-tbdEa?dl=0
% https://www.dropbox.com/sh/bm0at2dr81isk70/AABD67D_bF8A0NFa4gtt2dHNa?dl=0
% 
% Another highly helpful resource in this context is this Nature Methods
% paper by Tor and Wani Woo
% https://www.nature.com/articles/s41596-019-0289-5
%
% Lukas' version of the script for this paper can be found here
% https://www.dropbox.com/sh/v2nsgoqmbi0cqnk/AAD6I1Gn5KUM6aViom4TLeVJa?dl=0
% 
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   April, 2021
%__________________________________________________________________________
% @(#)% LaBGAS_run_MVPA_regression_single_trial     v1.0        
% last modified: 2021/04/01

%% LOAD FMRI_DATA_ST OBJECT AND (RE)DEFINE DIREECTORIES
%--------------------------------------------------------------------------
if ~exist('this_dat','var')
    load('single_trial_fmri_data_st_object.mat');
end

if ~exist('model','var')
model = 'model_2_CANlab_single_trial_GLM';
end

if ~exist('resultsdir','dir')
rootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS';
secondleveldir = fullfile(rootdir,'secondlevel\',model);
resultsdir = fullfile(secondleveldir,'results');
end


%% CLEAN DATA AND MASK
%--------------------------------------------------------------------------
% NOTE: contrary to @bogpetre in his scripts, I have already excluded the
% trials/runs/subjects without behavioral outcome in the previous script

% SET VIF THRESHOLD
vif_threshold = 5;

% DEFINE SUBJECT IDENTIFIERS PRIOR TO REMOVING BAD TRIALS
subject_id_vifs = this_dat.metadata_table.subject_id;
[uniq_subject_id_vifs, ~, subject_id_vifs] = unique(subject_id_vifs,'stable');
n_subj_vifs = size(uniq_subject_id_vifs,1);

% PLOT VIFS AND CALCULATE PERCENTAGE EXCEEDING THRESHOLD DEFINED ABOVE
% over subjects
v1=figure;
hold off
v1 = plot(this_dat.metadata_table.vif_values);
yline(5);
title('vif values over all trials');
xlabel('trial');
ylabel('variance inflation factor');

good_trials_idx = this_dat.metadata_table.vif_values < vif_threshold;
bad_trials_perc = sum(~good_trials_idx)./size(this_dat.metadata_table.vif_values,1).*100;
sprintf('%4.2f percent of trials exceeds a vif threshold of %d, indicating multicollinearity with noise regressors; script will remove them',bad_trials_perc,vif_threshold)

% per subject
v2=figure;
for i = 1:n_subj_vifs
    this_idx_vifs = find(i == subject_id_vifs);
    this_vifs = this_dat.metadata_table.vif_values(this_idx_vifs);
    
    this_good_trials_idx = this_vifs < vif_threshold;
    this_bad_trials_perc = sum(~this_good_trials_idx)./size(this_vifs,1).*100;
    sprintf('%4.2f percent of trials for subject %s exceeds a vif threshold of %d, indicating multicollinearity with noise regressors; script will remove them',this_bad_trials_perc,uniq_subject_id_vifs{i},vif_threshold)

    subplot(ceil(sqrt(n_subj_vifs)), ceil(n_subj_vifs/ceil(sqrt(n_subj_vifs))), i);
    hold off
    v2 = plot(this_vifs);
    yline(5);
    box off
    title(uniq_subject_id_vifs{i});
    xlabel('trial');
    ylabel('vif');
end

% REMOVE CON IMAGES CORRESPONDING TO TRIALS EXCEEDING VIF THRESHOLDS FROM
% DATA OBJECT
% NOTE: contrary to the standard CANlab fmri_data object, @bogpetre's
% fmri_data_st object applies the idx to all fields of the object,
% including the metadata_table etc, hence it should be used for all of our
% purposes, particularly single trial analysis
this_dat = this_dat.get_wh_image(good_trials_idx);

% MASK
gray_mask = fmri_mask_image('gray_matter_mask.img');
this_dat = this_dat.apply_mask(gray_mask);
% this_dat.mask = gray_mask; % fmri_data.apply_mask does not seem to update mask info of the object automatically, so we do that manually here
% this_dat.mask_descrip = 'default CANlab gray matter mask';

% ZSCORE BEHAVIORAL OUTCOME
% NOTE: useful for more interpretable values of prediction MSE
this_dat.Y = zscore(this_dat.Y);

% CLEAR SUBJECT IDENTIFIERS AND REDEFINE ON REDUCED DATA OBJECT
clear subject_id_vifs uniq_subject_id_vifs n_subj_vifs this_idx_vifs;

subject_id = this_dat.metadata_table.subject_id;
[uniq_subject_id, ~, subject_id] = unique(subject_id,'stable');
n_subj = size(uniq_subject_id,1);


%% DATA VISUALISATION PRIOR TO MODEL BUILDING
%--------------------------------------------------------------------------

% BETA IMAGES
h1=figure;
for j = 1:n_subj
    subj_idx = j == subject_id;
    this_subj_dat = this_dat.dat(:,subj_idx);
    q(j,:) = quantile(this_subj_dat(:),[0.025,0.5,0.975]);
    mu = mean(mean(this_subj_dat(:)));
    sd = std(this_subj_dat(:));
    h1 = plot([mu-sd, mu+sd],[j,j],'-');
    hold on;
    h2 = plot(mu,j,'o');
    h2.Color = h1.Color;
end
box off
title('Distribution of beta weights');
xlabel('\beta');
ylabel('Subject');
hold off

p = get(gcf,'Position');
set(gcf,'Position',[p(1:2),1024,2048]);

% BEHAVIORAL OUTCOME
% over subjects
b1=figure;
hold off;
b1=histogram(this_dat.Y);
box off
title('Histogram of single trial somatic symptom ratings');
xlabel('symptom rating');
ylabel('n(observations)');

% per subject
b2=figure;
for i = 1:n_subj
    this_idx_Y = find(i == subject_id);
    this_Y = this_dat.Y(this_idx_Y);

    subplot(ceil(sqrt(n_subj)), ceil(n_subj/ceil(sqrt(n_subj))), i);
    hold off
    b2 = histogram(this_Y);
    box off
    title(uniq_subject_id{i});
    xlabel('symptom rating');
    ylabel('n(obs)');
end


%% CROSS-VALIDATION FOLD SELECTION
%--------------------------------------------------------------------------
% NOTE: balancing over groups, stratifying over subjects (i.e. leave whole
% subject out)

group = this_dat.metadata_table.patient;
cv = cvpartition2(group, 'KFOLD', 5, 'Stratify', subject_id);
    fold_labels = zeros(size(this_dat.dat,2),1);
    for j = 1:cv.NumTestSets
        fold_labels(cv.test(j)) = j;
    end
    
    
%% FIT SINGLE-LEVEL MVPA MODELS ON RAW AND TRANSFORMED DATA
%--------------------------------------------------------------------------
% NOTE: we use support vector regression here, but can easily be changed to
% another algorithm for continuous outcomes included in CANlab's predict
% function, for example LASSO-PCR or PLS
% type help predict in Matlab command window for more info

this_dat_c = this_dat.rescale('centerimages');
this_dat_l2 = this_dat.rescale('l2norm_images');
this_dat_z = this_dat.rescale('zscoreimages');

% default
t0 = tic;
[d_cverr, d_stats, d_optout] = predict(this_dat, 'algorithm_name', 'cv_svr', ...
    'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
d_t = toc(t0);

fprintf('PCR r = %0.3f\n', corr(d_stats.yfit, this_dat.Y));

figure
line_plot_multisubject(this_dat.Y, d_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'SVR Estimated Symptoms','(cross validated)'})

figure
d_stats.weight_obj.montage;

% centered
t0 = tic;
[c_cverr, c_stats, c_optout] = predict(this_dat_c, 'algorithm_name', 'cv_svr', ...
    'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
c_t = toc(t0);

fprintf('PCR r = %0.3f\n', corr(c_stats.yfit, this_dat.Y));

figure
line_plot_multisubject(this_dat.Y, c_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'SVR Estimated Symptoms','(cross validated)'})

figure
c_stats.weight_obj.montage;

% l2normed
t0 = tic;
[l2_cverr, l2_stats, l2_optout] = predict(this_dat_l2, 'algorithm_name', 'cv_svr', ...
    'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
l2_t = toc(t0);

fprintf('PCR r = %0.3f\n', corr(l2_stats.yfit, this_dat.Y));

figure
line_plot_multisubject(this_dat.Y, l2_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'SVR Estimated Symptoms','(cross validated)'})

figure
l2_stats.weight_obj.montage;

% zscored
t0 = tic;
[z_cverr, z_stats, z_optout] = predict(this_dat_z, 'algorithm_name', 'cv_svr', ...
    'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
z_t = toc(t0);

fprintf('PCR r = %0.3f\n', corr(z_stats.yfit, this_dat.Y));

figure
line_plot_multisubject(this_dat.Y, z_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'SVR Estimated Symptoms','(cross validated)'})

figure
z_stats.weight_obj.montage;


%% FIT MULTILEVEL MVPA MODEL
%--------------------------------------------------------------------------

% DETERMINE MAXIMUM NUMBER OF COMPONENTS
max_comp = floor(size(this_dat.dat,2).*0.75 - n_subj);
% NOTE: we specify the maximum number of components as < the number of columns in
% this_dat.dat (n_subjects*n_conditions(in every fold)) to avoid overfitting in multilevel models, 
% where we need to leave df for the random intercepts (upper bound 1 df per random intercept hence subject)

% FIT SINGLE LEVEL MODEL
[pcr_cverr, pcr_stats,pcr_optout] = this_dat.predict('algorithm_name','cv_pcr',...
    'nfolds',fold_labels, 'numcomponents',max_comp);
fprintf('PCR r = %0.3f\n', corr(pcr_stats.yfit,this_dat.Y));

figure
line_plot_multisubject(this_dat.Y, pcr_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'PCR Estimated Symptoms','(cross validated)'})

figure
pcr_stats.weight_obj.montage;

% FIT MULTILEVEL MODEL W/FIXED PARAMETER ESTIMATION
% split maximum amount of components in between and within
n_bt_comp = floor(0.75*n_subj);
n_wi_comp = max_comp - n_bt_comp;
% NOTE: max between = n_subj IN EVERY FOLD (hence n_subj - 20% in 5-fold CV), 
% and you want to put more money on within since this typically explains
% more variance

% overall model prediction
[mlpcr_cverr, mlpcr_stats, mlpcr_optout] = this_dat.predict('algorithm_name','cv_mlpcr',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id);
fprintf('multilevel PCR r = %0.3f\n',corr(mlpcr_stats.yfit, this_dat.Y));
% lukasvo76: algorithm option 'cv_mlpcr' requires subject identifier,
% which makes sense since this is a multilevel/mixed model
% note that fold labels are the same, since they respect subject membership

figure
line_plot_multisubject(this_dat.Y, mlpcr_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'MLPCR Estimated Symptoms','(cross validated)'})

figure
mlpcr_stats.weight_obj.montage;

figure
subplot(1,2,1)
line_plot_multisubject(pcr_stats.yfit, mlpcr_stats.yfit, 'subjid', subject_id);
xlabel({'PCR model prediction'}); ylabel('Multilevel PCR model prediction');
axis square
subplot(1,2,2);
plot(pcr_optout{1}(:),mlpcr_optout{1}(:),'.');
lsline;
xlabel('PCR model weights'); ylabel('Multilevel PCR model weights');
axis square
% NOTE: contrary to @bogpetre's walkthrough, the
% pcr and the multilevel pcr models are not exactly equivalent anymore
% since I have been specifying the number of components

% get the variance explained by the between and within component
% NOTE: These functions call the same thing under the hood, 
% but simply perform cross validation using ONLY between or within
% subject models.
[mlpcr_bt_cverr, mlpcr_bt_stats] = this_dat.predict('algorithm_name','cv_mlpcr_bt',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
pred_bt = mlpcr_bt_stats.yfit;

[mlpcr_wi_cverr, mlpcr_wi_stats] = this_dat.predict('algorithm_name','cv_mlpcr_wi',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
pred_wi = mlpcr_wi_stats.yfit;
% lukasvo76: algorithm options are created by bogpetre

fprintf('Between subject PCR components r = %0.3f\n', corr(mlpcr_bt_stats.yfit, this_dat.Y));
fprintf('Within subject PCR components r = %0.3f\n', corr(mlpcr_wi_stats.yfit, this_dat.Y));

figure
subplot(1,2,1)
line_plot_multisubject(this_dat.Y, pred_bt, 'subjid', subject_id);
xlabel({'Observed pain'}); ylabel('Between subject components'' prediction');
axis square
subplot(1,2,2)
line_plot_multisubject(this_dat.Y, pred_wi, 'subjid', subject_id);
xlabel({'Observed pain'}); ylabel('Within subject components'' prediction');
axis square

% FIT MULITLEVEL MODEL W/ MIXED PARAMETER ESTIMATION
% NOTE: this is not part of the walkthrough yet, but @bogpetre pushed
% mlpcr3 function to CanlabCore
% bogpetre:
% main function is mlpcr3, so help mlpcr3 for usage options.
% basically it's the same as cv_mlpcr, except there's a randInt, randSlope and fitlmeOpts now
% fitlmeOpts get passed on to fitlme. It picks some sensible defaults
% randSlope makes things much slower. randInt is roughly the sae order of magnitude as running the fixed effects version
% the function defaults to fixed effects by default, so it's a drop in replacement
% for mlpcr2.m (aka cv_mlpcr)

% overall model prediction including random intercept only
[mlpcr3_cverr, mlpcr3_stats, mlpcr3_optout] = this_dat.predict('algorithm_name','cv_mlpcr3',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1);
fprintf('multilevel PCR r = %0.3f\n',corr(mlpcr3_stats.yfit, this_dat.Y));
% NOTE: compare with code in previous section and note change of algorithm_name option to cv_mlpcr3 and addition of randInt
% option - see help mlpcr3 for more details

figure
line_plot_multisubject(this_dat.Y, mlpcr3_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'MLPCR3 Estimated Symptoms','(cross validated)'})

figure
mlpcr3_stats.weight_obj.montage;

figure
subplot(1,2,1)
line_plot_multisubject(mlpcr_stats.yfit, mlpcr3_stats.yfit, 'subjid', subject_id);
xlabel({'Multilevel PCR model prediction'}); ylabel('Multilevel PCR model random int model prediction');
axis square
subplot(1,2,2);
plot(mlpcr_optout{1}(:),mlpcr3_optout{1}(:),'.');
lsline;
xlabel('Multilevel PCR model weights'); ylabel('Multilevel PCR model random int model weights');
axis square
% NOTE: note that models are very similar but not exactly equivalent

% overall model prediction including random intercept and random slope
[mlpcr3rs_cverr, mlpcr3rs_stats, mlpcr3rs_optout] = this_dat.predict('algorithm_name','cv_mlpcr3',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1, 'randSlope', 1);
fprintf('multilevel PCR r = %0.3f\n',corr(mlpcr3rs_stats.yfit, this_dat.Y));

figure
line_plot_multisubject(this_dat.Y, mlpcr3rs_stats.yfit, 'subjid', subject_id);
xlabel({'Observed Symptoms','(average over conditions)'}); ylabel({'MLPCR3rs Estimated Symptoms','(cross validated)'})

figure
mlpcr3rs_stats.weight_obj.montage;


%% SAVE STATS FOR ALL MODELS
%--------------------------------------------------------------------------

savefilename = fullfile(resultsdir, 'single_trial_MVPA_results.mat');
save(savefilename, 'cv','d_stats','pcr_stats','mlpcr_stats','mlpcr3_stats','mlpcr3rs_stats', '-v7.3');