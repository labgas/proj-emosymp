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
% LOAD DATA

% conditions
idx = ~isnan(DAT.BEHAVIOR.behavioral_data_table.NA_neg);
behdat = DAT.BEHAVIOR.behavioral_data_table(idx,:);
behdat = sortrows(behdat,'patient','ascend');
behdat.patient(behdat.patient == -1) = 2;
group = [behdat.patient; behdat.patient; behdat.patient];
condition = [ones(height(behdat),1); 2.*ones(height(behdat),1); 3.*ones(height(behdat),1)];
NA = [behdat.NA_neg; behdat.NA_neu; behdat.NA_pos];
symptoms = [behdat.symptoms_neg; behdat.symptoms_neu; behdat.symptoms_pos];
outcomes = {NA,symptoms};
outcome_names = {'negative affect','symptoms'};

for o = 1:size(outcomes,2)
    D{o} = [outcomes{o},condition,group];
    for i = 1:3
        for j = 1:2
        data{i,j} = D{o}(D{o}(:, 2) == i & D{o}(:, 3) ==j);
        end
    end
    
    data_all{o} = data;
end

% contrasts
for k = 1:max(unique(behdat.patient))
    
    NA_contrast{k} = [behdat.NA_neg_neu(behdat.patient==k),behdat.NA_neg_pos(behdat.patient==k)];
    symptoms_contrast{k} = [behdat.symptoms_neg_neu(behdat.patient==k),behdat.symptoms_pos_neu(behdat.patient==k)];
    
end

outcomes_contrast = {NA_contrast,symptoms_contrast};
contrast_names = {'negative versus neutral','negative versus positive'};


% CONDITIONS - RM RAINCLOUD PLOT
% individual figures
for o = 1:size(outcomes,2)
    f  = figure('Position', fig_position,'WindowState','maximized');
    h   = rm_raincloud(data_all{o}, cl, 0, 'rash');
    ax{o} = gca;
    ax{o}.FontSize = 14;
    ax{o}.FontName = 'Cambria';
    ax{o}.XAxis.LineWidth = 1;
    ax{o}.YAxis.LineWidth = 1;
    % title(['Figure 1' newline 'Repeated measures raincloud plot']);
    xlabel({strcat(outcome_names{o},' rating'),''},'FontSize',24,'FontWeight','bold');
    ylabel({'','condition'},'FontSize',24,'FontWeight','bold');
    yticklabels({'\fontsize{20} \bf positive','\fontsize{20} \bf neutral','\fontsize{20} \bf negative'});
    legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
    for i = 1:3
        for j = 1:2
            h.s{i, j}.SizeData = 150;
        end
    end
    
    % save
    print(f,fullfile(figspubdir,strcat('behav_',outcome_names{o},'.png')),'-dpng','-r600');
    
    f_all{o}= f;
    h_all{o} = h;
    
    clear f h;
end

% integrated two-panel figure for publication
fa  = figure('WindowState','maximized');
subplot(1,2,1)
ha   = rm_raincloud(data_all{1}, cl, 0, 'rash');
ax1 = gca;
ax1.FontSize = 12;
ax1.FontName = 'Cambria';
ax1.XAxis.LineWidth = 1;
ax1.YAxis.LineWidth = 1;
ax1.Position = [0.07 0.1100 0.42 0.8150];
% title(['Figure 1' newline 'Repeated measures raincloud plot']);
xlabel({strcat(outcome_names{1},' rating'),''},'FontSize',20,'FontWeight','bold');
ylabel({'','condition'},'FontSize',20,'FontWeight','bold');
yticklabels({'\fontsize{16} \bf positive','\fontsize{16} \bf neutral','\fontsize{16} \bf negative'});
% legend([h.l(1,1) h.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold');
    for i = 1:3
        for j = 1:2
            ha.s{i, j}.SizeData = 100;
        end
    end
subplot(1,2,2)
hb   = rm_raincloud(data_all{2}, cl, 0, 'rash');
ax2 = gca;
ax2.FontSize = 12;
ax2.FontName = 'Cambria';
ax2.XAxis.LineWidth = 1;
ax2.YAxis.LineWidth = 1;
ax2.Position = [0.5703 0.1100 0.42 0.8150];
% title(['Figure 1' newline 'Repeated measures raincloud plot']);
xlabel({strcat(outcome_names{2},' rating'),''},'FontSize',20,'FontWeight','bold');
ylabel({'','condition'},'FontSize',20,'FontWeight','bold');
yticklabels({'\fontsize{16} \bf positive','\fontsize{16} \bf neutral','\fontsize{16} \bf negative'});
legend([hb.l(1,1) hb.l(1,2)],{'FSS patients','healthy controls'},'Location','best','FontSize',20,'FontWeight','bold','Box','off');
    for i = 1:3
        for j = 1:2
            hb.s{i, j}.SizeData = 100;
        end
    end
    
% fa.Position = [-150 50 2250 971]; % can be used to change position of
% entire figure
    
print(fa,fullfile(figspubdir,strcat('behav_',outcome_names{1},'_',outcome_names{2},'.png')),'-dpng','-r600');


% CONTRASTS - CLASSIC RAINCLOUD PLOT

for o = 1:size(outcomes,2)
    
    for m = 1:size(NA_contrast,2)
        
        f2 = figure('Position', fig_position,'WindowState','maximized');
        h1 = raincloud_plot(outcomes_contrast{o}{1}(:,m), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 1, 'line_width', 3);
        h2 = raincloud_plot(outcomes_contrast{o}{2}(:,m), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
        h1{1}.EdgeColor = 'none';
        h2{1}.EdgeColor = 'none';
        h1{2}.SizeData = 50;
        h2{2}.SizeData = 50;
        ax3{o,m} = gca;
        ax3{o,m}.FontSize = 14;
        ax3{o,m}.FontName = 'Cambria';
        ax3{o,m}.XAxis.LineWidth = 1;
        ax3{o,m}.YAxis.LineWidth = 1;
        xlabel({'',[outcome_names{o},' rating']},'FontSize',24,'FontWeight','bold');
        legend([h1{1} h2{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
        title(contrast_names{m},'FontSize',28,'FontWeight','bold');
        ylim([(h2{3}.Position(2)+0.10.*h2{3}.Position(2)) (max([h1{1}.YData h2{1}.YData])+0.05.*max([h1{1}.YData h2{1}.YData]))]);
        box off
    
    
        % save
        print(f2,fullfile(figspubdir,strcat('behav_contrasts_',outcome_names{o},'_',contrast_names{m},'.png')),'-dpng','-r600');

        f2_all{o,m}= f2;
        h1_all{o,m} = h1;
        h2_all{o,m} = h2;

        clear f2 h1 h2;
        
    end
    
end


%% SIGNATURE RESPONSES
%--------------------------------------------------------------------------
% LOAD DATA
sigdat = DAT.SIG_contrasts.raw.dotproduct;
behdat_full = DAT.BEHAVIOR.behavioral_data_table;
behdat_full.patient(behdat_full.patient == -1) = 2;

for n = 1:max(unique(behdat_full.patient))
    
    NPS{n} = [sigdat.NPS.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPS.Negative_v_Positive(behdat_full.patient==n)];
    NPSpos{n} = [sigdat.NPSpos.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPSpos.Negative_v_Positive(behdat_full.patient==n)];
    NPSneg{n} = [sigdat.NPSneg.Negative_v_Neutral(behdat_full.patient==n),sigdat.NPSneg.Negative_v_Positive(behdat_full.patient==n)];
    PINES{n} = [sigdat.PINES.Negative_v_Neutral(behdat_full.patient==n),sigdat.PINES.Negative_v_Positive(behdat_full.patient==n)];
    
end

signatures = {NPS, NPSpos, NPSneg, PINES};
signature_names = {'NPS', 'NPS positive', 'NPS negative', 'PINES'};
contrast_names = {'negative versus neutral','negative versus positive'};

npsposregions = DAT.NPSsubregions.npspos_by_region_contrasts;
npsposregions_names = DAT.NPSsubregions.posnames;
npsposregions_names_full = {'vermis','R insula','R V1','R thalamus','L insula','R dorsal posterior insula','R S2 operculum','anterior midcingulate'};

for t = 1:max(unique(behdat_full.patient))
    
    for u = 1:size(npsposregions{1},2)
       regions{u,t} = [npsposregions{1}(behdat_full.patient==t,u),npsposregions{2}(behdat_full.patient==t,u)];
    end
    
end

% INDIVIDUAL FIGURES
for s = 1:size(signatures,2)
    
    for p = 1:size(NPS,2)
        
        f3 = figure('Position', fig_position,'WindowState','maximized');
        h3 = raincloud_plot(signatures{s}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
             'box_col_match', 1, 'line_width', 3);
        h4 = raincloud_plot(signatures{s}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 3);
        h3{1}.EdgeColor = 'none';
        h4{1}.EdgeColor = 'none';
        h3{2}.SizeData = 50;
        h4{2}.SizeData = 50;
        ax4{s,p} = gca;
        ax4{s,p}.FontSize = 14;
        ax4{s,p}.FontName = 'Cambria';
        ax4{s,p}.XAxis.LineWidth = 1;
        ax4{s,p}.YAxis.LineWidth = 1;
        xlabel({'',[signature_names{s},' response']},'FontSize',24,'FontWeight','bold');
        legend([h3{1} h4{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',24,'FontWeight','bold','Box','off');
        title(contrast_names{p},'FontSize',28,'FontWeight','bold');
        ylim([(h4{3}.Position(2)+0.10.*h4{3}.Position(2)) (max([h3{1}.YData h4{1}.YData])+0.05.*max([h3{1}.YData h4{1}.YData]))]);
        box off
    
    
        % save
        print(f3,fullfile(figspubdir,strcat(signature_names{s},'_',contrast_names{p},'.png')),'-dpng','-r600');

        f3_all{s,p}= f3;
        h3_all{s,p} = h3;
        h4_all{s,p} = h4;

        clear f3 h3 h4;
        
    end
    
end

clear p s;


% INTEGRATED PANEL FIGURES FOR PUBLICATION

% PINES and NPS

f4  = figure('Position', fig_position,'WindowState','maximized');

for p = 1:size(NPS,2)
    subplot(2,2,p)
        h5 = raincloud_plot(signatures{4}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h6 = raincloud_plot(signatures{4}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h5{1}.EdgeColor = 'none';
        h6{1}.EdgeColor = 'none';
        h5{2}.SizeData = 50;
        h6{2}.SizeData = 50;
        ax4{p} = gca;
        ax4{p}.FontSize = 12;
        ax4{p}.FontName = 'Cambria';
        ax4{p}.XAxis.LineWidth = 1;
        ax4{p}.YAxis.LineWidth = 1;
        ax4{p}.YAxisLocation = 'origin';
        ax4{p}.YTick = [];
        ax4{p}.LineWidth = 0.25;
        xlabel({[signature_names{4},' response']},'FontSize',18,'FontWeight','bold');
        legend([h5{1} h6{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
        title({(contrast_names{p}),''},'FontSize',20,'FontWeight','bold');
        ylim([(h6{3}.Position(2)+0.10.*h6{3}.Position(2)) (max([h5{1}.YData h6{1}.YData])+0.05.*max([h5{1}.YData h6{1}.YData]))]);
        box off
end

for q = 3:(size(NPS,2)+2)
    subplot(2,2,q)
        h7 = raincloud_plot(signatures{1}{1}(:,q-2), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h8 = raincloud_plot(signatures{1}{2}(:,q-2), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h7{1}.EdgeColor = 'none';
        h8{1}.EdgeColor = 'none';
        h7{2}.SizeData = 50;
        h8{2}.SizeData = 50;
        ax4{q} = gca;
        ax4{q}.FontSize = 12;
        ax4{q}.FontName = 'Cambria';
        ax4{q}.XAxis.LineWidth = 1;
        ax4{q}.YAxis.LineWidth = 1;
        ax4{q}.YAxisLocation = 'origin';
        ax4{q}.YTick = [];
        ax4{q}.LineWidth = 0.25;
        xlabel({[signature_names{1},' response']},'FontSize',18,'FontWeight','bold');
        legend([h7{1} h8{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
%         title(contrast_names{q-2},'FontSize',20,'FontWeight','bold');
        ylim([(h8{3}.Position(2)+0.10.*h8{3}.Position(2)) (max([h7{1}.YData h8{1}.YData])+0.05.*max([h7{1}.YData h8{1}.YData]))]);
        box off
end
    
print(f4,fullfile(figspubdir,strcat(signature_names{4},'_',signature_names{1},'.png')),'-dpng','-r600');

clear p q;

% NPS positive and negative

f5  = figure('Position', fig_position,'WindowState','maximized');

for p = 1:size(NPS,2)
    subplot(2,2,p)
        h9 = raincloud_plot(signatures{2}{1}(:,p), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h10 = raincloud_plot(signatures{2}{2}(:,p), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h9{1}.EdgeColor = 'none';
        h10{1}.EdgeColor = 'none';
        h9{2}.SizeData = 50;
        h10{2}.SizeData = 50;
        ax4{p} = gca;
        ax4{p}.FontSize = 12;
        ax4{p}.FontName = 'Cambria';
        ax4{p}.XAxis.LineWidth = 1;
        ax4{p}.YAxis.LineWidth = 1;
        ax4{p}.YAxisLocation = 'origin';
        ax4{p}.YTick = [];
        ax4{p}.LineWidth = 0.25;
        xlabel({[signature_names{2},' response']},'FontSize',18,'FontWeight','bold');
        legend([h9{1} h10{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
        title({(contrast_names{p}),''},'FontSize',20,'FontWeight','bold');
        ylim([(h10{3}.Position(2)+0.10.*h10{3}.Position(2)) (max([h9{1}.YData h10{1}.YData])+0.05.*max([h9{1}.YData h10{1}.YData]))]);
        box off
end

for q = 3:(size(NPS,2)+2)
    subplot(2,2,q)
        h11 = raincloud_plot(signatures{3}{1}(:,q-2), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                     'box_col_match', 1);
        h12 = raincloud_plot(signatures{3}{2}(:,q-2), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
             'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1);
        h11{1}.EdgeColor = 'none';
        h12{1}.EdgeColor = 'none';
        h11{2}.SizeData = 50;
        h12{2}.SizeData = 50;
        ax4{q} = gca;
        ax4{q}.FontSize = 12;
        ax4{q}.FontName = 'Cambria';
        ax4{q}.XAxis.LineWidth = 1;
        ax4{q}.YAxis.LineWidth = 1;
        ax4{q}.YAxisLocation = 'origin';
        ax4{q}.YTick = [];
        ax4{q}.LineWidth = 0.25;
        xlabel({[signature_names{3},' response']},'FontSize',18,'FontWeight','bold');
        legend([h11{1} h12{1}], {'FSS patients', 'healthy controls'},'Location','best','FontSize',14,'FontWeight','bold','Box','off');
%         title(contrast_names{q-2},'FontSize',20,'FontWeight','bold');
        ylim([(h12{3}.Position(2)+0.10.*h12{3}.Position(2)) (max([h11{1}.YData h12{1}.YData])+0.05.*max([h11{1}.YData h12{1}.YData]))]);
        box off
end
    
print(f5,fullfile(figspubdir,strcat(signature_names{2},'_',signature_names{3},'.png')),'-dpng','-r600');

clear p q;


% NPS subregions

for c = 1:size(contrast_names,2)
    
    f{c}  = figure('Position', fig_position,'WindowState','maximized');
    title(contrast_names{c},'FontSize',24,'FontWeight','bold');

    for p = 1:size(regions,1)
       subplot(3,3,p)
            h13 = raincloud_plot(regions{p,1}(:,c), 'box_on', 1, 'color', cb(4,:), 'alpha', 0.5,...
                         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
                         'box_col_match', 1, 'line_width', 1.5);
            h14 = raincloud_plot(regions{p,2}(:,c), 'box_on', 1, 'color', cb(1,:), 'alpha', 0.5,...
                 'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 1, 'line_width', 1.5);
            h13{1}.EdgeColor = 'none';
            h14{1}.EdgeColor = 'none';
            h13{2}.SizeData = 20;
            h14{2}.SizeData = 20;
            ax5{p} = gca;
            ax5{p}.FontSize = 12;
            ax5{p}.FontName = 'Cambria';
            ax5{p}.XAxis.LineWidth = 1;
            ax5{p}.YAxis.LineWidth = 1;
            ax5{p}.YAxisLocation = 'origin';
            ax5{p}.YTick = [];
            ax5{p}.LineWidth = 0.25;
            xlabel({[npsposregions_names_full{p},' response']},'FontSize',14,'FontWeight','bold');
            legend([h13{1} h14{1}], {'FSS', 'controls'},'Location','best','FontSize',10,'FontWeight','bold','Box','off');
    %         title({(npsposregions_names{p}),''},'FontSize',20,'FontWeight','bold');
            ylim([(h14{3}.Position(2)+0.10.*h14{3}.Position(2)) (max([h13{1}.YData h14{1}.YData])+0.05.*max([h13{1}.YData h14{1}.YData]))]);
            box off
    end
    
    sgtitle({contrast_names{c},''},'FontSize',20,'FontWeight','bold','FontName','Cambria');
    
    print(f{c},fullfile(figspubdir,strcat('npssubregions_',contrast_names{c},'.png')),'-dpng','-r600');
    
end