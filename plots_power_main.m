close all;clear all;clc;
%% get data
PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

windows={'base_beta','stim_beta','prespeech_beta','speech_beta','rebound_beta'}; 
windows_low={'base_lowbeta','stim_lowbeta','prespeech_lowbeta','speech_lowbeta','rebound_lowbeta'};
windows_high={'base_highbeta','stim_highbeta','prespeech_highbeta','speech_highbeta','rebound_highbeta'};
ii=1:length(SUBJECTS);
xpos=nan(1,numel(windows));
xpos1=nan(1,numel(windows));
xpos2=nan(1,numel(windows));
xpos3=nan(1,numel(windows));
timing_windows=struct();


%% beta low high
fig1=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

for c=1:2
if c==1;subplot(2,1,1);else;subplot(2,1,2);end
for w=1:numel(windows)
    window=string(windows(w));
    window_low=string(windows_low(w));
    window_high=string(windows_high(w));
    disp(strcat('ANALYSIS :'," ",window))  
    %subplot(2,numel(windows),w)
    
    
    power_ecog=[];powerlow_ecog=[];powerhigh_ecog=[];
    power_dbs=[];powerlow_dbs=[];powerhigh_dbs=[];
    for i=ii
        if i==6;continue;end
        %if i==14;continue;end
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

        PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');

        % take sessions with DBS LFP data
        session=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
        id_session=session.id(strcmpi(session.type, 'LEAD'));
    
        % upload useful annot tables
        coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding')); coding=coding(coding.session_id==id_session,:);
        cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise')); cue=cue(cue.session_id==id_session, :);

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

        
        switch window
            case 'base_beta'
                window_start= mean(tf_baseline_t0(:,1),'omitnan');
                window_end= mean(tf_baseline_t0(:,2),'omitnan');
                timing_windows.baseline.start=window_start;
                timing_windows.baseline.end=window_end;
            case 'stim_beta'
                window_start= x_stim1;
                window_end= x_stim3_end;
                timing_windows.stimulus.start=window_start;
                timing_windows.stimulus.end=window_end;
            case 'prespeech_beta'
                window_start= x_stim3_end;
                window_end= x_syl1;
                timing_windows.prespeech.start=window_start;
                timing_windows.prespeech.end=window_end;
            case 'speech_beta'
                window_start= x_syl1;
                window_end= x_syl3_end;
                timing_windows.speech.start=window_start;
                timing_windows.speech.end=window_end;
            case 'rebound_beta'
                window_start= tf_rebound_t0(1);
                window_end= tf_rebound_t0(2);
                timing_windows.rebound.start=window_start;
                timing_windows.rebound.end=window_end;
        end


        tab=readtable(strcat('annot/general CTAR/',SUBJECT,'_speech_pow_db.txt')); 
        power_ecog=[power_ecog ; tab.(window)(startsWith(tab.label,'ecog'))]; % array
        powerlow_ecog=[powerlow_ecog ; tab.(window_low)(startsWith(tab.label,'ecog'))]; % array
        powerhigh_ecog=[powerhigh_ecog ; tab.(window_high)(startsWith(tab.label,'ecog'))]; % array
        power_dbs=[power_dbs ; tab.(window)(startsWith(tab.label,'dbs_L'),:)]; % array
        powerlow_dbs=[powerlow_dbs ; tab.(window_low)(startsWith(tab.label,'dbs_L'),:)]; % array
        powerhigh_dbs=[powerhigh_dbs ; tab.(window_high)(startsWith(tab.label,'dbs_L'),:)]; % array

    end
    xpos(w)=window_end-(window_end-window_start)/2;
    xpos1(w)=xpos(w)-0.2;
    xpos2(w)=xpos(w);
    xpos3(w)=xpos(w)+0.2;

    if c==1
        powEcog.(window)=power_ecog;
        powEcog.(window_low)=powerlow_ecog;
        powEcog.(window_high)=powerhigh_ecog;

        boxplot(powEcog.(window_low),'Positions',xpos1(w));ylim([-6,4]);hold on
        boxplot(powEcog.(window_high),'Positions',xpos2(w));ylim([-6,4]);
        boxplot(powEcog.(window),'Positions',xpos3(w));ylim([-6,4]);
        value=mean(powEcog.(window),'omitnan');
        text(xpos3(w)+0.1, value, sprintf('%.3f', value), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        ylimits=ylim;
        marker_pos = ylimits(1) + 0.8*(ylimits(2)-ylimits(1)); 
        [h, p, ci, stats] = ttest2(powEcog.(window_low), powEcog.(window_high));
    else
        powDbs.(window)=power_dbs;
        powDbs.(window_low)=powerlow_dbs;
        powDbs.(window_high)=powerhigh_dbs;

        boxplot(powDbs.(window_low),'Positions',xpos1(w));ylim([-6,4]);hold on
        boxplot(powDbs.(window_high),'Positions',xpos2(w));ylim([-6,4]);
        boxplot(powDbs.(window),'Positions',xpos3(w));ylim([-6,4]);
        value=mean(powDbs.(window),'omitnan');
        text(xpos3(w)+0.1, value, sprintf('%.3f', value), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        ylimits=ylim;
        marker_pos = ylimits(1) + 0.8*(ylimits(2)-ylimits(1)); 
        [h, p, ci, stats] = ttest2(powDbs.(window_low), powDbs.(window_high));
    end
    % plot a marker if significant comparison is found
    if p < 0.05
        if p < 0.001;n_markers=3;
        elseif p < 0.01;n_markers=2;
        else; n_markers=1;
        end
        plot([xpos1(w),xpos2(w)], [marker_pos, marker_pos], 'k-')
        markers=repmat('*',1,n_markers);
        text(mean([xpos1(w),xpos2(w)]), marker_pos + 0.2 ,markers,'HorizontalAlignment', 'center', 'FontSize',12)
    end

end

% colors
h = findobj(gca,'Tag','Box'); 
for box=1:length(h) % start from the end to count boxes!
    if mod(box,3)==0
        patch(get(h(box),'XData'),get(h(box),'YData'),[0.5 1 0.5], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    elseif mod(box,3)==2
        patch(get(h(box),'XData'),get(h(box),'YData'),[0 0.5 0], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    else
        patch(get(h(box), 'XData'), get(h(box), 'YData'), [1 0 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    end
end

ylim([-6,4]);
ylimits=ylim;
ylabel('Power [dB]')
fill([timing_windows.baseline.start timing_windows.baseline.end timing_windows.baseline.end timing_windows.baseline.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [1 0.8 0.4], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([timing_windows.stimulus.start timing_windows.stimulus.end timing_windows.stimulus.end timing_windows.stimulus.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [0.3 0.4 0.9], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([timing_windows.prespeech.start timing_windows.prespeech.end timing_windows.prespeech.end timing_windows.prespeech.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)],[0.5 0.8 0.4], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([timing_windows.speech.start timing_windows.speech.end timing_windows.speech.end timing_windows.speech.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [1 0.3 0.3], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([timing_windows.rebound.start timing_windows.rebound.end timing_windows.rebound.end timing_windows.rebound.start], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], [1 0.6 0.4], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
set(gca, 'XTick', xpos, 'XTickLabel', {'baseline','stimulus','prespeech','speech','postspeech'});

if c==1;title('ECoG');else;title('DBS');end
end
sgtitle({'Low vs high beta power per time windows','',''})

%%
fig2=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
tit2='Power per beta band levels and window';
sgtitle({tit2,' '});
for p=1:(2*numel(windows))
    subplot(2,numel(windows),p)
    if p<=numel(windows)
        data_ecog=[powEcog.(windows{p}), powEcog.(windows_low{p}), powEcog.(windows_high{p})];
        boxplot(data_ecog, 'Labels',{'Beta general', 'Low beta', 'High beta'})
        ylim([-15,5])
        title(strrep(windows{p},'_',' '))

    else
        data_dbs=[powDbs.(windows{p-numel(windows)}), powDbs.(windows_low{p-numel(windows)}), powDbs.(windows_high{p-numel(windows)})];
        boxplot(data_dbs, 'Labels',{'Beta general', 'Low beta', 'High beta'})
        ylim([-15,5])
        title(strrep(windows{p-numel(windows)},'_',' '))
    end
end


% NON PRENDE SGTITLE
clear power_ecog powEcog

%% ecog location

locations={'PMC','POC','IFC','SMC','TPOJ','LTC','DLPFC','IFOC','AAC','IPC'};
for w=1:numel(windows)
    window=string(windows(w));
    disp(strcat('ANALYSIS :'," ",window))  

    for l=1:numel(locations)
        location=locations{l};
        power_ecog.(location)=[];
        for i=ii
            if i==6;continue;end
            if i==14;continue;end
            SUBJECT=strcat('DBS',string(SUBJECTS(i)));
            disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
            % paths 
            PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
            electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
            cfg=[];
            cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                        % threshold_subj or threshold_weight
            electrode=bml_getEcogArea(cfg,electrode);
            
            channels=electrode.electrode(strcmp(electrode.HCPMMP1_area,location));

            tab=readtable(strcat('annot/general CTAR/',SUBJECT,'_speech_pow_db.txt')); 
            power_ecog.(location)=[power_ecog.(location) ; tab.(window)(find(ismember(tab.label,channels)))]; % array    
        end
        powEcog.(window).(location)=power_ecog.(location);
    end
end

fig3=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
for w=1:numel(windows)
    window=string(windows(w));
    subplot(1,numel(windows),w)
    for l=1:numel(locations)
        location=locations{l};
        boxplot(powEcog.(window).(location),'Positions',l);
        ylim([-15,5])
        hold on;
    end
    set(gca, 'XTick', 1:numel(locations), 'XTickLabel', locations);
    title(strrep(window,'_',' '))
end
tit3='Power per ECoG location';
sgtitle({tit3,' '})

%%
disp('---saving')
saveas(fig1,strcat('images/general CTAR/cortex/','power low high','.png'))
saveas(fig1,strcat('images/general CTAR/cortex/','power low high','.fig'))

saveas(fig2,strcat('images/general CTAR/cortex/',tit2,'.png'))
saveas(fig2,strcat('images/general CTAR/cortex/',tit2,'.fig'))

saveas(fig3,strcat('images/general CTAR/cortex/',tit3,'.png'))
saveas(fig3,strcat('images/general CTAR/cortex/',tit3,'.fig'))

%%
close all



