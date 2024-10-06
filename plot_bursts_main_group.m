close all;clear all;clc;
%%
PATH_DATA='Z:\DBS';

DATE=datestr(now,'yyyymmdd');
format long
% from previous analysis:
A=2; % median distance stim1 onset - syl 1 onset
B=1.26; % median distance syl1 onset - syl3 offset
IT=1.93; % median inter trial distance (syl3 offset - stim 1 onset)
LEFT=1.5;
RIGHT=1.5; 

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

%% 
ii=1:nuumel(SUBJECTS);

for i=ii
    %open a subject dir
    if i==6;continue;end
    if i==14;continue;end
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');

    % take sessions with DBS LFP data
    session=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    id_session=session.id(strcmpi(session.type, 'LEAD'));

    % upload useful annot tables
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding')); coding=coding(coding.session_id==id_session,:);
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise')); cue=cue(cue.session_id==id_session, :);
    prod_triplet=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_produced_triplet'));prod_triplet=prod_triplet(prod_triplet.session_id==id_session,:);
    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    electrode=electrode(:, {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1','HCPMMP1_weight_1'});

    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);

    tokeepCod=(~isnan(coding.syl1_onset));
    tokeepCue=(~isnan(cue.stim1_starts));
    switch string(SUBJECTS(i))
        case {'3003', '3010', '3015', '3020', '3025', '3027'};tokeepCue=tokeepCod;
        case {'3012','3014'};tokeepCod=tokeepCue;
        case {'3024'};tokeepCod=tokeepCod & tokeepCue;tokeepCue=tokeepCod;
        case {'3021'};idx=find(~tokeepCue)+1;idx=idx(1:2);tokeepCod(idx)=0;
    end
    center=coding.syl1_onset(tokeepCod);

    x=center - cue.stim1_starts(tokeepCue); %TRIAL SPECIFIC
    tf_baseline_t0=[-(x+0.9) -(x+0.1)];

    x=mean(coding.syl3_offset(tokeepCod) - center,'omitnan'); %MEDIATING
    tf_rebound_t0=[x+0.1 x+0.9];

    x_stim1=median(cue.stim1_starts(tokeepCue)-center,'omitnan');
    x_stim3_end=median(cue.stim3_ends(tokeepCue)-center,'omitnan');
    x_syl1=median(coding.syl1_onset(tokeepCod)-center,'omitnan');
    x_syl3_end=median(coding.syl3_offset(tokeepCod)-center,'omitnan');
    
    windows={'baseline','stimulus','prespeech','speech','rebound'};
    timing_windows=struct();
    for w=1:numel(windows)
        window=windows{w};
        switch window
            case 'baseline'
                timing_windows.(window).start= mean(tf_baseline_t0(:,1),'omitnan');
                timing_windows.(window).end= mean(tf_baseline_t0(:,2),'omitnan');
            case 'stimulus'
                timing_windows.(window).start= x_stim1;
                timing_windows.(window).end= x_stim3_end;
            case 'prespeech'
                timing_windows.(window).start= x_stim3_end;
                timing_windows.(window).end= x_syl1;
            case 'speech'
                timing_windows.(window).start= x_syl1;
                timing_windows.(window).end= x_syl3_end;
            case 'rebound'
                timing_windows.(window).start= tf_rebound_t0(1);
                timing_windows.(window).end= tf_rebound_t0(2);
        end
    end

    tab_beta=readtable(strcat('annot/general CTAR/',SUBJECT,'_beta_bursts.txt'));
    tab_fooof=readtable(strcat('annot/general CTAR/',SUBJECT,'_fooof_bursts.txt'));
    tab_cycle=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt'));

    % filter out bursts <=100ms NB CAMBIA MOLTO SE SI VEDE SOLO SOTTO 100ms!!
    tab_beta=tab_beta(tab_beta.duration>0.1,:);
    tab_fooof=tab_fooof(tab_fooof.duration>0.1,:);
    tab_cycle=tab_cycle(tab_cycle.duration>0.1,:);

    electrodes_subject= unique([tab_beta.chan_id ; tab_fooof.chan_id ; tab_cycle.chan_id]);
    electrodes_areas= unique( electrode.HCPMMP1_area (~cellfun('isempty',electrode.HCPMMP1_area)) );
    electrodes_areas{end+1} = 'dbs_L';

    %% comparison with 
    %% power measures comparison
    continue
    disp('--start power measurements analysis')

    pType={'power','pow_max','pow_mean'};
    vType={'band_amp','volt_peak','volt_amp'};
   
    for a=1:numel(electrodes_areas)
        area= electrodes_areas{a};
        fig1=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

        if strcmp(area,'dbs_L')
            sgtitle({strcat(SUBJECT," ",num2str(sum(startsWith(unique(tab_beta.chan_id),area)))," ", area), ' '})
            tab_beta_el=tab_beta( startsWith(tab_beta.chan_id,area),:);
            tab_fooof_el=tab_fooof( startsWith(tab_fooof.chan_id,area),:);
            tab_cycle_el=tab_cycle( startsWith(tab_cycle.chan_id,area),:);
        else
            area_name= string( electrode.HCPMMP1_area_fullname( strcmp(electrode.HCPMMP1_area,area)) );
            sgtitle({strcat(SUBJECT," ",num2str(numel(area_name))," ",area_name(1)), ' '})
            tab_beta_el=tab_beta( ismember( tab_beta.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            tab_fooof_el=tab_fooof( ismember( tab_fooof.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            tab_cycle_el=tab_cycle( ismember( tab_cycle.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
        end

        for p=1:numel(pType)
            xpower=pType{p};
            if p==3;ypower=pType{1};else;ypower=pType{p+1};end
            xvoltage=vType{p};
            if p==3;yvoltage=vType{1};else;yvoltage=vType{p+1};end
        
            subplot(3,3,p)
            %[~, d] = ksdensity([tab_beta.(xpower),tab_beta.(ypower)]);
            scatter(tab_beta_el.(xpower),tab_beta_el.(ypower),'filled', 'MarkerFaceAlpha', 0.2)
            xlabel(strrep(xpower,'_',' ')),ylabel(strrep(ypower,'_',' '))
            title('Beta general')
        
            subplot(3,3,3+p)
            scatter(tab_fooof_el.(xpower),tab_fooof_el.(ypower),'filled', 'MarkerFaceAlpha', 0.2)
            xlabel(strrep(xpower,'_',' ')),ylabel(strrep(ypower,'_',' '))
            title('Fooof')
        
            subplot(3,3,6+p)
            scatter(tab_cycle_el.(xvoltage),tab_cycle_el.(yvoltage),'filled', 'MarkerFaceAlpha', 0.2)
            xlabel(strrep(xvoltage,'_',' ')),ylabel(strrep(yvoltage,'_',' '))
            title('Bycycle')
        end
        
        saveas(fig1,strcat('images/areas/',SUBJECT," ",area," ",'power measures.png'));
        saveas(fig1,strcat('images/areas/',SUBJECT," ",area," ",'power measures.fig'));
        close all

    end
      
    
    %% micellaneous comparison
    disp('--start micellaneous analysis')
    % I keep pow_max and volt_amp

    for a=1:numel(electrodes_areas)
        area= electrodes_areas{a};
        fig2=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

        if strcmp(area,'dbs_L')
            sgtitle({strcat(SUBJECT," ",num2str(sum(startsWith(unique(tab_beta.chan_id),area)))," ",area), ' '})
            tab_beta_el=tab_beta( startsWith(tab_beta.chan_id,area),:);
            tab_fooof_el=tab_fooof( startsWith(tab_fooof.chan_id,area),:);
            tab_cycle_el=tab_cycle( startsWith(tab_cycle.chan_id,area),:);
        else
            area_name= string( electrode.HCPMMP1_area_fullname( strcmp(electrode.HCPMMP1_area,area)) );
            sgtitle({strcat(SUBJECT," ",num2str(numel(area_name))," ",area_name(1)), ' '})
            tab_beta_el=tab_beta( ismember( tab_beta.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            tab_fooof_el=tab_fooof( ismember( tab_fooof.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            tab_cycle_el=tab_cycle( ismember( tab_cycle.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
        end
    
        subplot(2,3,1)
        scatter(tab_beta_el.('pow_max'),tab_beta_el.('duration'),'filled', 'MarkerFaceAlpha', 0.2)
        xlabel('max power'),ylabel('duration')
        title({'Correlations between max power and duration','Beta general'})
        
        subplot(2,3,2)
        scatter(tab_fooof_el.('pow_max'),tab_fooof_el.('duration'),'filled', 'MarkerFaceAlpha', 0.2)
        xlabel('max power'),ylabel('duration')
        title({'Correlations between max power and duration','Fooof'})
        
        subplot(2,3,3)
        scatter(tab_cycle_el.('volt_amp'),tab_cycle_el.('duration'),'filled', 'MarkerFaceAlpha', 0.2)
        xlabel('band voltage'),ylabel('duration')
        title({'Correlations between max power and duration','Bycycle'})
    
        % cycles - duration comparison
        subplot(2,3,4)
        scatter(tab_cycle_el.('cycles'),tab_cycle_el.('duration'),'filled', 'MarkerFaceAlpha', 0.2)
        xlabel('cycles'),ylabel('duration')
        title({'Comparison cycles - duration','Bycycle'})
        
        % fooof freq vs cycles 
        subplot(2,3,5)
        histogram(tab_cycle_el.frequency,'FaceColor',[0.4,0.6,1]);hold on
        iFreq=mean(tab_fooof_el.ifreq);
        xline([iFreq iFreq],'LineWidth',2,'Color',[1,0.7,0.4])
        xrange = [iFreq - 2.5, iFreq + 2.5];
        x_fill = [xrange(1), xrange(1), xrange(2), xrange(2)];
        y_fill = [min(ylim), max(ylim), max(ylim), min(ylim)];
        fill(x_fill, y_fill, [1,0.7,0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        title('Bycycle frequecies vs fooof iFreq')

        saveas(fig2,strcat('images/areas/',SUBJECT," ",area," ",'micellaneous.png'));
        saveas(fig2,strcat('images/areas/',SUBJECT," ",area," ",'micellaneous.fig'));
        close all
    end
    
    %% overlap general
    % uso IoU, se tanti zeri considera GIoU
    
    time = linspace(-(A+LEFT), B+RIGHT, 500 );

    for a=1:numel(electrodes_areas)
        disp('--start Jaccard Index analysis')
        area= electrodes_areas{a};
        fig3=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

        if strcmp(area,'dbs_L')
            sgtitle({strcat(SUBJECT," ",num2str(sum(startsWith(unique(tab_beta.chan_id),area)))," ",area), ' '})
            tab_beta_el=tab_beta( startsWith(tab_beta.chan_id,area),:);
            tab_fooof_el=tab_fooof( startsWith(tab_fooof.chan_id,area),:);
            tab_cycle_el=tab_cycle( startsWith(tab_cycle.chan_id,area),:);
        else
            area_name= string( electrode.HCPMMP1_area_fullname( strcmp(electrode.HCPMMP1_area,area)) );
            sgtitle({strcat(SUBJECT," ",num2str(numel(area_name))," ",area_name(1)), ' '})
            tab_beta_el=tab_beta( ismember( tab_beta.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            tab_fooof_el=tab_fooof( ismember( tab_fooof.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            tab_cycle_el=tab_cycle( ismember( tab_cycle.chan_id, electrode.electrode(strcmp(electrode.HCPMMP1_area,area))),:);
            electrodes_analysis= intersect( electrodes_subject, electrode.electrode(strcmp(electrode.HCPMMP1_area,area)));
        end

        tabgt=tab_beta_el;

        trial_ids=unique(tabgt.trial_id);
        for m=1:2
            subplot(4,3,m)
            switch m
                case 1;tabcmp=tab_cycle_el;tit='Jaccard Index per windows - bycycle';
                case 2;tabcmp=tab_fooof_el;tit='Jaccard Index per windows - fooof';
            end
            for w=1:numel(windows)
                % mask the window
                start_window= timing_windows.(windows{w}).start;
                end_window= timing_windows.(windows{w}).end;
                [~,idx_start]=min(abs(time-start_window));
                [~,idx_end]=min(abs(time-end_window));
                window_time=time(idx_start:idx_end);
                iou=nan(1);
                for e=1:numel(electrodes_analysis)
                    % mark for each trial
                    iou_trials=nan(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        trial_id=trial_ids(j);
                        subtabgt=tabgt( strcmp(tabgt.chan_id,electrodes_analysis(e)) & tabgt.trial_id==trial_id  & strcmp(tabgt.window,windows(w)),:);
                        subtabcmp=tabcmp( strcmp(tabcmp.chan_id,electrodes_analysis(e)) & tabcmp.trial_id==trial_id & strcmp(tabcmp.window,windows(w)) ,:);
                        
                        maskgt=false(size(window_time));
                        maskcmp=false(size(window_time));
                        try
                            for pp=1:height(subtabgt)
                                [~,idx_start]=min(abs(window_time-subtabgt.starts(pp)));
                                [~,idx_end]=min(abs(window_time-subtabgt.ends(pp)));
                                maskgt(idx_start:idx_end)=1;
                            end     
                            for pp=1:height(subtabcmp)
                                [~,idx_start]=min(abs(window_time-subtabcmp.starts(pp)));
                                [~,idx_end]=min(abs(window_time-subtabcmp.ends(pp)));
                                maskcmp(idx_start:idx_end)=1;
                            end
                        catch
                        end
                        iou_trials(j) = sum( maskgt & maskcmp ) / sum( maskgt | maskcmp );
                    end
                    iou=[iou , iou_trials];
                end
                boxplot(iou,'Positions',w)
                ylim([0,1])
                hold on
            end
            set(gca, 'XTick', 1:numel(windows), 'XTickLabel', windows);
            title(tit)
        end

        subplot(4,3,3)
        for m=1:2
            switch m
                case 1;tabcmp=tab_cycle_el;
                case 2;tabcmp=tab_fooof_el;
            end
            iou=nan(1);
            for e=1:numel(electrodes_analysis)
                iou_trials=nan(1,numel(trial_ids));
                for j=1:numel(trial_ids)
                    trial_id=trial_ids(j);
                    subtabgt=tabgt( strcmp(tabgt.chan_id,electrodes_analysis(e)) & tabgt.trial_id==trial_id ,:);
                    subtabcmp=tabcmp( strcmp(tabcmp.chan_id,electrodes_analysis(e)) & tabcmp.trial_id==trial_id ,:);
                    
                    maskgt=false(size(time));
                    maskcmp=false(size(time));
                    try
                        for pp=1:height(subtabgt)
                            [~,idx_start]=min(abs(time-subtabgt.starts(pp)));
                            [~,idx_end]=min(abs(time-subtabgt.ends(pp)));
                            maskgt(idx_start:idx_end)=1;
                        end 
                        for pp=1:height(subtabcmp)
                            [~,idx_start]=min(abs(time-subtabcmp.starts(pp)));
                            [~,idx_end]=min(abs(time-subtabcmp.ends(pp)));
                            maskcmp(idx_start:idx_end)=1;
                        end
                    catch
                    end
                    iou_trials(j) = sum( maskgt & maskcmp ) / sum( maskgt | maskcmp );
                end
                iou=[iou , iou_trials];
            end
            boxplot(iou,'Positions',m)
            ylim([0,1])
            hold on
        end
        title('Jaccard Index overall')
        set(gca, 'XTick', 1:2, 'XTickLabel', {'bycycle','fooof'});

        % burst probability
        disp('--start burst probability analysis')
        for m=1:3
            subplot(4,3,3+m)
            switch m
                case 1;tab=tab_beta_el;tit='Burst probability per window - beta general';
                case 2;tab=tab_fooof_el;tit='Burst probability per window - fooof';
                case 3;tab=tab_cycle_el;tit='Burst probability per window - bycycle';  
            end
            for w=1:numel(windows)
                % mask the window
                start_window= timing_windows.(windows{w}).start;
                end_window= timing_windows.(windows{w}).end;
                [~,idx_start]=min(abs(time-start_window));
                [~,idx_end]=min(abs(time-end_window));
                window_time=time(idx_start:idx_end);
                prob=nan(1);
                for e=1:numel(electrodes_analysis)
                    % mask for each trial
                    prob_trials=nan(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        trial_id=trial_ids(j);
                        subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_id  & strcmp(tab.window,windows(w)),:);
                        
                        mask=false(size(window_time));
                        try
                            for pp=1:height(subtab)
                                [~,idx_start]=min(abs(window_time-subtab.starts(pp)));
                                [~,idx_end]=min(abs(window_time-subtab.ends(pp)));
                                mask(idx_start:idx_end)=1;
                            end 
                        catch
                        end
                        prob_trials(j) = sum(mask) / length(mask);
                    end
                    prob = [prob , prob_trials];
                end
                boxplot(prob,'Positions',w)
                ylim([0,1])
                hold on
            end
            set(gca, 'XTick', 1:numel(windows), 'XTickLabel', windows);
            title(tit)
        end   
    
        % duration
        toEvaluate='duration';
        disp('--start burst duration analysis')
        for m=1:3
            subplot(4,3,6+m)
            switch m
                case 1;tab=tab_beta_el;tit='Burst duration per window - beta general';
                case 2;tab=tab_fooof_el;tit='Burst duration per window - fooof';
                case 3;tab=tab_cycle_el;tit='Burst duration per window - bycycle';  
            end
            for w=1:numel(windows)
                % mask the window
                start_window= timing_windows.(windows{w}).start;
                end_window= timing_windows.(windows{w}).end;
                [~,idx_start]=min(abs(time-start_window));
                [~,idx_end]=min(abs(time-end_window));
                window_time=time(idx_start:idx_end);
                dur=nan(1);
                for e=1:numel(electrodes_analysis)
                    % costruisco una maschera per ogni trial
                    dur_trials=cell(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        trial_id=trial_ids(j);
                        subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_id  & strcmp(tab.window,windows(w)),:);
                        measure=nan(1,height(subtab));
                        try
                            for pp=1:height(subtab)
                                measure(pp)=subtab.(toEvaluate);
                            end 
                        catch
                        end
                        dur_trials{j} = measure(pp);
                    end
                    dur = [dur , cell2mat(dur_trials)];
                end
                boxplot(dur,'Positions',w)
%                 ylim([0,1])
                hold on
            end
            set(gca, 'XTick', 1:numel(windows), 'XTickLabel', windows);
            title(tit)
        end 

        % power
        disp('--start burst power analysis')
        for m=1:3
            subplot(4,3,9+m)
            switch m
                case 1;tab=tab_beta_el;tit='Burst power per window - beta general';toEvaluate='pow_max';
                case 2;tab=tab_fooof_el;tit='Burst power per window - fooof';toEvaluate='pow_max';
                case 3;tab=tab_cycle_el;tit='Burst voltage per window - bycycle';toEvaluate='volt_amp';
            end
            for w=1:numel(windows)
                % mask the window
                start_window= timing_windows.(windows{w}).start;
                end_window= timing_windows.(windows{w}).end;
                [~,idx_start]=min(abs(time-start_window));
                [~,idx_end]=min(abs(time-end_window));
                window_time=time(idx_start:idx_end);

                pow=nan(1);
                for e=1:numel(electrodes_analysis)
                    % mask for each trial
                    pow_trials=cell(1,numel(trial_ids));
                    for j=1:numel(trial_ids)
                        trial_id=trial_ids(j);
                        subtab=tab( strcmp(tab.chan_id,electrodes_analysis(e)) & tab.trial_id==trial_id  & strcmp(tab.window,windows(w)),:);
                        measure=nan(1,height(subtab));
                        try
                            for pp=1:height(subtab)
                                measure(pp)=subtab.(toEvaluate);
                            end 
                        catch
                        end
                        pow_trials{j} = measure(pp);
                    end        
                    pow=[pow , cell2mat(pow_trials)];
                end
                boxplot(pow,'Positions',w)
                %ylim([0,1])
                hold on
            end
            set(gca, 'XTick', 1:numel(windows), 'XTickLabel', windows);
            title(tit)
        end
        
        saveas(fig3,strcat('images/areas/',SUBJECT," ",area," ",'pow_volt.png'));
        saveas(fig3,strcat('images/areas/',SUBJECT," ",area," ",'pow_volt.fig'));
        close all
    
    end
end  
disp('finish :)')


