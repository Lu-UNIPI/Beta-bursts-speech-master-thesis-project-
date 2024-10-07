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
measures={'duration','volt_amp','time_rdsym','time_ptsym','frequency','n_bursts','perc_bursts'}; % n burst and probability
areas=unique(tab_areas.area);
windows={'baseline','stimulus','prespeech','speech','rebound'};


%% get data 
ii=14%7:numel(SUBJECTS);

%     windows={'baseline','rebound'}; 
%     feature_ecog_baseline = [];
%     feature_ecog_rebound = [];
%     feature_dbs_baseline = [];
%     feature_dbs_rebound = [];
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
    annot_subject=table(); row=1;
        
    % get BURSTS DATA
    tab=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt'));
    tab_occurrence=readtable(strcat('annot/general CTAR/bursts occurrence/',SUBJECT,'_bursts_occurrence.txt'));
    electrodes_analysis=unique(tab.chan_id);

    for e=1:numel(electrodes_analysis)
        if startsWith(electrodes_analysis{e},'dbs_R');continue;end
        tab_electrode=tab(strcmp(tab.chan_id,electrodes_analysis{e}),:);
        tab_occurrence_electrode=tab_occurrence(strcmp(tab_occurrence.label,electrodes_analysis{e}),:);
        if i==14;area=electrode.HCPMMP1_area(strcmp(electrode.electrode(1:130,:),electrodes_analysis(e)));
            else;area=electrode.HCPMMP1_area(strcmp(electrode.electrode,electrodes_analysis(e)));
        end
        if (isempty(area) && startsWith(electrodes_analysis{e},'dbs_L'));area={'dbs'};end

        for m=1:numel(measures)
            measure=measures{m};
 
            annot_subject.id(row)=row;
            annot_subject.label(row)=electrodes_analysis(e);
            annot_subject.area(row)=area;
            annot_subject.measure(row)=measures(m);
            
            switch measure
                case {'duration','volt_amp','time_rdsym','time_ptsym','frequency'}
                feature_baseline= tab_electrode.(measure)(strcmp(tab_electrode.window50,'baseline')); % array
                case {'n_bursts','perc_bursts'}
                feature_baseline= tab_occurrence_electrode.(measure)(strcmp(tab_occurrence_electrode.window50,'baseline')); % array
            end

            for w=2:numel(windows)
                window=windows{w};
                switch measure
                    case {'duration','volt_amp','time_rdsym','time_ptsym','frequency'}
                    feature_comparison= tab_electrode.(measure)(strcmp(tab_electrode.window50,window)); % array
                    case {'n_bursts','perc_bursts'}
                    feature_comparison= tab_occurrence_electrode.(measure)(strcmp(tab_occurrence_electrode.window50,window)); % array
                end

                if ~isempty(feature_baseline) && ~isempty(feature_comparison)
                    warning off
                    try 
                        [p,h,stats] = ranksum(feature_comparison, feature_baseline);
                        if isfield(stats, 'zval')
                            zval=stats.zval;
                        else
                            n1=numel(feature_comparison);
                            n2=numel(feature_baseline);
                            R1=stats.ranksum;
                            mu_R = n1 * (n1 + n2 + 1) / 2;
                            sigma_R = sqrt(n1 * n2 * (n1 + n2 + 1) / 12);
                            zval = (R1 - mu_R) / sigma_R;
                        end
                        annot_subject.(strcat('wtest_z_',window))(row)=zval;
                        annot_subject.(strcat('wtest_pv_',window))(row)=p;
                    catch
                    end
                    
    
                    [~,p,~,stats] = ttest2(feature_comparison, feature_baseline);
                    annot_subject.(strcat('ttest_t_',window))(row)=stats.tstat;
                    annot_subject.(strcat('ttest_pv_',window))(row)=p;
                    warning on
 
                end
            end

            row=row+1;
        end
    end
    writetable(annot_subject, strcat('annot/general CTAR/areas/',SUBJECT," ",'bycycle features comparison.txt'), 'Delimiter', '\t');
    disp('-----------------------------------')
end






