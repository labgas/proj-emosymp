%% DEFINE PATHS AND VARIABLES TO BE USED IN ALL ANALYSES
% -------------------------------------------------------------------------
addpath(genpath('C:\Users\lukas\Documents\GitHub\proj-emosymp')); % add proj-emosymp Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\RainCloudPlots')); % add RainCloudPlots Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\GitHub\Robust_Statistical_Toolbox')); % add Robust Statistical Toolbox Github repo to your path
addpath(genpath('C:\Users\lukas\Documents\MATLAB\cbrewer')); % add colorbrewer to your path for more color options

a_emosymp_m1_s1_set_up_paths_always_run_first

figspubdir = fullfile(resultsdir,'figures_publication');
if ~isfolder(figspubdir)
    mkdir(figspubdir);
end

load(fullfile(resultsdir,'image_names_and_setup.mat'));

try
    % get nice colours from colorbrewer
    % (https://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
catch
    % if you don't have colorbrewer, accept these far more boring colours
    cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
end

cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);

fig_position = [200 200 600 400]; % coordinates for figures


%% BEHAVIORAL DATA
%--------------------------------------------------------------------------
idx = ~isnan(DAT.BEHAVIOR.behavioral_data_table.NA_neg);
behdat = DAT.BEHAVIOR.behavioral_data_table(idx,:);
behdat = sortrows(behdat,'patient','ascend');
behdat.patient(behdat.patient == -1) = 2;
group = [behdat.patient; behdat.patient; behdat.patient];
condition = [ones(height(behdat),1); 2.*ones(height(behdat),1); 3.*ones(height(behdat),1)];
NA = [behdat.NA_neg; behdat.NA_neu; behdat.NA_pos];
symptoms = [behdat.symptoms_neg; behdat.symptoms_neu; behdat.symptoms_pos];
outcomes = {NA,symptoms};
outcome_names = {'NA','symptoms'};

for o = 1:size(outcomes,2)
    D{o} = [outcomes{o},condition,group];
    for i = 1:3
        for j = 1:2
        data{i,j} = D{o}(D{o}(:, 2) == i & D{o}(:, 3) ==j);
        end
    end
    
    data_all{o} = data;
end

for o = 1:size(outcomes,2)
    f  = figure('Position', fig_position,'WindowState','maximized');
    h   = rm_raincloud(data_all{o}, cl, 0, 'rash');
    % title(['Figure 1' newline 'Repeated measures raincloud plot']);
    xlabel({strcat(outcome_names{o},' rating'),''},'FontSize',24,'FontWeight','bold');
    ylabel({'','condition'},'FontSize',24,'FontWeight','bold');
    yticklabels({'\fontsize{20} positive','\fontsize{20} neutral','\fontsize{20} negative'});
    legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','northeast','FontSize',24,'FontWeight','bold');

    for i = 1:3
        for j = 1:2
            h.s{i, j}.SizeData = 150;
        end
    end
    
    % save
    print(f,fullfile(figspubdir,strcat('behav_',outcome_names{o},'.png')),'-dpng');
    
    f_all{o}= f;
    h_all{o} = h;
    
    clear f h;
end

for k = 1:max(unique(behdat.patient))
    
    NA_contrast{k} = [behdat.NA_neg_neu(behdat.patient==k),behdat.NA_neg_pos(behdat.patient==k)];
    symptoms_contrast{k} = [behdat.symptoms_neg_neu(behdat.patient==k),behdat.symptoms_pos_neu(behdat.patient==k)];
    
end

outcomes_contrast = {NA_contrast,symptoms_contrast};
contrast_names = {'negative versus neutral','negative versus positive'};

for o = 1:size(outcomes,2)
    
    for m = 1:size(NA_contrast,2)
        
        f2 = figure('Position', fig_position,'WindowState','maximized');
        h1 = raincloud_plot(outcomes_contrast{o}{1}(:,m), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 0);
        h2 = raincloud_plot(outcomes_contrast{o}{2}(:,m), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
        h1{2}.SizeData = 50;
        h2{2}.SizeData = 50;
        xlabel({'',[contrast_names{m},' ',outcome_names{o},' rating']},'FontSize',16,'FontWeight','bold');
        legend([h1{1} h2{1}], {'FSS patients', 'healthy controls'},'Location','northeast','FontSize',16,'FontWeight','bold');
    %     title(['Figure M7' newline 'A) Dodge Options Example 1']);
        set(gca,'YLim', [-0.6 1.5]);
        box off
    
    
        % save
        print(f2,fullfile(figspubdir,strcat('behav_contrasts_',outcome_names{o},'_',contrast_names{m},'.png')),'-dpng');

        f2_all{o,m}= f2;
        h1_all{o,m} = h1;
        h2_all{o,m} = h2;

        clear f2 h1 h2;
        
    end
    
end