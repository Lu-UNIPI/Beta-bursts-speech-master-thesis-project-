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
measures={'duration','volt_amp','time_rdsym','time_ptsym','frequency'}; 
%measures={'n_bursts','perc_bursts'};
areas=unique(tab_areas.area); %{'PMC','SMC','IPC'};



for m=1:numel(measures)
for a=1:numel(areas)
    % -----------------PUT HERE WHAT TO PLOT!!-------------------------

    measure=measures{m};
    toEvaluate=measure;

    area=areas{a};

    comparison='rebound'; % compare vs baseline

    % -----------------------------------------------------------------

    %% get data 
    ii=1:numel(SUBJECTS);

%     windows={'baseline','rebound'}; 
%     feature_ecog_baseline = [];
%     feature_ecog_rebound = [];
%     feature_dbs_baseline = [];
%     feature_dbs_rebound = [];
    for i=ii
        %open a subject dir
        if i==6;continue;end
        if i==14;continue;end
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
        tab=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt')); 
        %tab=readtable(strcat('annot/general CTAR/bursts occurrence/',SUBJECT,'_bursts_occurrence.txt')); % cambia 3 label con chan_id
        
        area_channels=electrode.electrode(strcmp(electrode.HCPMMP1_area,area));
        tab_area=tab( ismember(tab.chan_id,area_channels) ,:);

        fig1=figure('units','normalized','outerposition',[0.03 0.03 0.95 0.95]);
        sgtitle(strcat(SUBJECT," ",'Comparison in bycycle features between baseline and'," ",comparison," - ",strrep(measure,'_',' ')," - ",area),'FontWeight', 'bold')
        electrodes_analysis=unique(tab_area.chan_id);
        numelect=numel(electrodes_analysis);
        if numelect<1;continue
        elseif numelect<12;hor=2;lung=ceil(numelect/2);
        else;hor=3;lung=ceil(numelect/3);
        end

        for e=1:numel(electrodes_analysis)
            tab_electrode=tab_area(strcmp(tab_area.chan_id,electrodes_analysis{e}),:);

            feature_baseline= tab_electrode.(toEvaluate)(strcmp(tab_electrode.window50,'baseline')); % array
            feature_comparison= tab_electrode.(toEvaluate)(strcmp(tab_electrode.window50,comparison)); % array

            subplot(hor,lung,e)
            histogram(feature_baseline,'FaceColor','g','FaceAlpha',0.4);hold on
            if strcmp(comparison,'rebound');col_comparison='r';
            elseif strcmp(comparison,'speech');col_comparison=[0.1,0.4,1];
            end
            histogram(feature_comparison,'FaceColor',col_comparison,'FaceAlpha',0.4)
            legend('baseline',comparison);
            xlabel(strrep(measure,'_',' '))
            title(strrep(electrodes_analysis{e},'_',' '))
            
            if any(~(isempty(feature_baseline)), ~(isempty(feature_comparison)))
                [p, h,stats] = ranksum(feature_comparison, feature_baseline);
                x_lim=xlim; y_lim=ylim;      
                if (p<0.05 && h==1);color='yellow';else;color='white';end
                text(x_lim(1) + 0.7*(x_lim(2)-x_lim(1)), y_lim(1) + 0.7*(y_lim(2)-y_lim(1)), sprintf('ranksum test \np= %.3f',p),...
                    'FontSize',10, 'FontWeight','bold','BackgroundColor', color, 'EdgeColor', 'black')
            end
            
            
        end
        saveas(fig1,strcat('images/bursts CTAR/areas/',area," ",SUBJECT," ",'baseline-',comparison," ",measure,'.png'))
        saveas(fig1,strcat('images/bursts CTAR/areas/',area," ",SUBJECT," ",'baseline-',comparison," ",measure,'.fig'))
        close all


    end

    disp('-----------------------------------')
    
end
end



