close all;clear all;clc;
%%
PATH_DATA='Z:\DBS';

DATE=datestr(now,'yyyymmdd');
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

tab_areas=readtable("HCPMMP1toAreas.txt");

%% analisi
% measures={'probability','duration','volt_amp','time_rdsym','time_ptsym','power'}; 
measures={'duration','volt_amp','time_rdsym','time_ptsym','frequency','n_bursts','perc_bursts'};
areas=unique(tab_areas.area); %{'PMC','SMC','IPC'};
areas(11)={'dbs'};

for m=1:numel(measures)
measure=measures{m};
disp('-----------------------------------')
disp(strcat('-- ',upper(measure)))

% -----------------PUT HERE WHAT TO PLOT!!-------------------------
comparison='speech'; % compare vs baseline
stat_toEvaluate='ttest_t_speech';
pv_toEvaluate='ttest_pv_speech';
% -----------------------------------------------------------------

fig1=figure('units','normalized','outerposition',[0.03 0.03 0.95 0.95]);
sgtitle(strcat('Comparison in bycycle features between baseline and'," ",comparison," - ",upper(strrep(measure,'_',' '))),'FontWeight', 'bold')

p_area=0;
for a=1:numel(areas)   
    area=areas{a};
    disp(strcat('-- ',area))
    if strcmp(area,'LTC') ||strcmp(area,'IFOC');continue;end
    p_area=p_area+1;
    subplot(3,3,p_area)

    %% get data 
    ii=1:numel(SUBJECTS);
    for i=ii
        %open a subject dir
        if i==6;continue;end
        %if i==14;continue;end
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
        
        PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
        electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
        electrode=electrode(:, {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1','HCPMMP1_weight_1'});
        cfg=[];
        cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                    % threshold_subj or threshold_weight
        electrode=bml_getEcogArea(cfg,electrode);
            
        % get BURSTS DATA
        tab_stats=readtable(strcat('annot/general CTAR/areas/',SUBJECT," ",'bycycle features comparison.txt'));
        
        if strcmp(area,'dbs')
            tab_stats_area=tab_stats( (startsWith(tab_stats.label,'dbs') & strcmp(tab_stats.measure,measure)) ,:);
            switch measure
                case {'duration','volt_amp','time_rdsym','time_ptsym','frequency'}
                    tab_meas=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt'));
                    tab_meas_area=tab_meas( (startsWith(tab_meas.chan_id,'dbs')) ,:);
                case {'n_bursts','perc_bursts'}
                    tab_meas=readtable(strcat('annot/general CTAR/bursts occurrence/',SUBJECT,'_bursts_occurrence.txt'));
                    tab_meas_area=tab_meas( (startsWith(tab_meas.label,'dbs')) ,:);
            end
        else
            area_channels=electrode.electrode(strcmp(electrode.HCPMMP1_area,area)); 
            tab_stats_area=tab_stats( (ismember(tab_stats.label,area_channels) & strcmp(tab_stats.measure,measure)) ,:);
            switch measure
                case {'probability','duration','volt_amp','time_rdsym','time_ptsym','frequency'}
                    tab_meas=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt'));
                    tab_meas_area=tab_meas( (ismember(tab_meas.chan_id,area_channels)) ,:);
                case {'n_bursts','perc_bursts'}
                    tab_meas=readtable(strcat('annot/general CTAR/bursts occurrence/',SUBJECT,'_bursts_occurrence.txt'));
                    tab_meas_area=tab_meas( (ismember(tab_meas.label,area_channels)) ,:);
            end
        end
        
        electrodes_analysis=unique(tab_stats_area.label);
        numelect=numel(electrodes_analysis);
        all_feature_baseline = nan;
        all_feature_comparison = nan;
        for e=1:numel(electrodes_analysis)
            if startsWith(electrodes_analysis,'dbs_R');continue;end
            try 
                tab_meas_electrode=tab_meas_area(strcmp(tab_meas_area.label,electrodes_analysis{e}),:);
            catch
                tab_meas_electrode=tab_meas_area(strcmp(tab_meas_area.chan_id,electrodes_analysis{e}),:);    
            end
            tab_stats_electrode=tab_stats_area(strcmp(tab_stats_area.label,electrodes_analysis{e}),:);

            feature_baseline= mean( tab_meas_electrode.(measure)(strcmp(tab_meas_electrode.window50,'baseline')) ,'omitnan'); % point
            feature_comparison= mean( tab_meas_electrode.(measure)(strcmp(tab_meas_electrode.window50,comparison)), 'omitnan'); % point
            stat_comparison= tab_stats_electrode.(stat_toEvaluate);
            pv_comparison= tab_stats_electrode.(pv_toEvaluate);
            
            if stat_comparison>0 
                if pv_comparison<0.001;color=[0.6, 0, 0];size=70; % dark red
                elseif pv_comparison<0.01;color=[1, 0.6, 0.6];size=50; % red
                elseif pv_comparison<0.05;color=[1, 0.8, 0.8];size=30; % light pink
                else; color=[0.9,0.9,0.9];size=10;
                end
            else
                if pv_comparison<0.001;color=[0, 0, 0.6];size=70; % dark blue
                elseif pv_comparison<0.01;color=[0.3, 0.6, 1];size=50; % blue
                elseif pv_comparison<0.05;color=[0.7, 0.9, 1];size=30; % light blue
                else; color=[0.9,0.9,0.9];size=10;
                end
            end     

            scatter(feature_baseline, feature_comparison,size,color,'filled', 'MarkerEdgeColor', 'none')
            hold on
%             xlim([min_value max_value]);
%             ylim([min_value max_value]);
            all_feature_baseline = [all_feature_baseline, feature_baseline];
            all_feature_comparison = [all_feature_comparison, feature_comparison];
        end
        min_val = min([nanmin(all_feature_baseline), nanmin(all_feature_comparison)]);
        max_val = max([nanmax(all_feature_baseline), nanmax(all_feature_comparison)]);
        plot([min_val max_val], [min_val max_val], 'k');
        axis square;
        xlabel('baseline')
        ylabel(comparison)
        title(area)
    end    
end
QUI
saveas(fig1,strcat('images/bursts CTAR/areas/group/','Baseline-',stat_toEvaluate," ",measure,'.png'))
saveas(fig1,strcat('images/bursts CTAR/areas/group/','Baseline-',stat_toEvaluate," ",measure,'.fig'))
close all
disp('-----------------------------------')
disp('-----------------------------------')
disp('-----------------------------------')

end



