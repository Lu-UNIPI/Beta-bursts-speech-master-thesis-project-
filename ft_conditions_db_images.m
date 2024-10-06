close all;clear all;clc;
%% 
%%% this script get raw data and 
% 1) open cleaned and epoched data


PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long
% from previous analysis:
A=2; % median distance stim1 onset - syl 1 onset
B=1.26; % median distance syl1 onset - syl3 offset
IT=1.93; % median inter trial distance (syl3 offset - stim 1 onset)
LEFT=1.5;
RIGHT=1.5; 
TF_RATE = 20; %Hz Sampling rate of time frequency plot
TF_FOI = round(10.^(0.30:0.05:2.4),2,'signif');
TF_BASELINE_WIDTH = 0.5;
beta=[11,22,36];

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for cc=1
for i=5%1:numel(SUBJECTS)  
    if i==6;continue;end
    %open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    % paths
    NAME_RAW_SIGNAL=strcat(SUBJECT,'_ft_raw_session');
    PATH_SIGNAL=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\FieldTrip\');
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');

    % take sessions with DBS LFP data
    session=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    id_session=session.id(strcmpi(session.type, 'LEAD'));

    % upload useful annot tables
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding')); coding=coding(coding.session_id==id_session,:);
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise')); cue=cue(cue.session_id==id_session, :);
    prod_triplet=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_produced_triplet'));
    prod_triplet=prod_triplet(prod_triplet.session_id==id_session,:);
    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    electrode=electrode(:, {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1','HCPMMP1_weight_1'});
    
    tokeepCod=(~isnan(coding.syl1_onset));
    tokeepCue=(~isnan(cue.stim1_starts));
    switch string(SUBJECTS(i))
        case {'3003', '3010', '3015', '3020', '3025', '3027'};tokeepCue=tokeepCod;
        case {'3014','3012'};tokeepCod=tokeepCue;
        case {'3024'};tokeepCod=tokeepCod & tokeepCue;tokeepCue=tokeepCod;
        case {'3021'};idx=find(~tokeepCue)+1;idx=idx(1:2);tokeepCod(idx)=0;
    end
    
    % load data
    load(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_syl.mat"));a=-(A+LEFT);b=B+RIGHT;center=coding.syl1_onset(tokeepCod);
    %load(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_stim.mat"));a=-(IT+LEFT);b=A+RIGHT;center=cue.stim1_starts(tokeepCue);
    
    % define epoch, baseline and rebound for ft
    tf_epoch_t0=E.time{:};tf_epoch_t0=[tf_epoch_t0(1) tf_epoch_t0(end)];
    tf_epoch_t0=round(tf_epoch_t0 .* TF_RATE) ./ TF_RATE;

    x=center - cue.stim1_starts(tokeepCue); %TRIAL SPECIFIC
    % x=median(coding.syl1_onset - cue.stim1_starts,'omitnan'); %MEDIATING
    tf_baseline_t0=[-(x+0.9) -(x+0.1)];
    tf_baseline_t0=round(tf_baseline_t0 .* TF_RATE) ./ TF_RATE;

    % x=coding.syl3_offset(tokeepCod) - center; %TRIAL SPECIFIC
    x=mean(coding.syl3_offset(tokeepCod) - center,'omitnan'); %MEDIATING
    tf_rebound_t0=[x+0.1 x+0.9];
    tf_rebound_t0=round(tf_rebound_t0 .* TF_RATE) ./ TF_RATE;

    % define stim and syl x
    x_stim1=median(cue.stim1_starts(tokeepCue)-center,'omitnan');
    x_stim2=median(cue.stim2_starts(tokeepCue)-center,'omitnan');
    x_stim3=median(cue.stim3_starts(tokeepCue)-center,'omitnan');
    x_stim3_end=median(cue.stim3_ends(tokeepCue)-center,'omitnan');
    x_stim1_end=median(cue.stim1_ends(tokeepCue)-center,'omitnan');
    x_stim2_end=median(cue.stim2_ends(tokeepCue)-center,'omitnan');
    x_syl1=median(coding.syl1_onset(tokeepCod)-center,'omitnan');
    x_syl2=median(coding.syl2_onset(tokeepCod)-center,'omitnan');
    x_syl3=median(coding.syl3_onset(tokeepCod)-center,'omitnan');
    x_syl1_end=median(coding.syl1_offset(tokeepCod)-center,'omitnan');
    x_syl2_end=median(coding.syl2_offset(tokeepCod)-center,'omitnan');
    x_syl3_end=median(coding.syl3_offset(tokeepCod)-center,'omitnan');

    %% select ecog and referencing
    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);
    
    % ecog
    choice='ecog';
    switch choice
        case 'area'
            channels=electrode.HCPMMP1_area;
            channels=channels(~(cellfun('isempty',channels)));
            idx=find(cellfun(@(x) strcmp(x,'AAC'),channels));
            % "PMC","POC","IFC","SMC","TPOJ","LTC","DLPFC","IFOC","AAC","IPC"
            if sum(idx)==0
                warning('No channels in this area');
            else
                cfg=[];
                cfg.channel     =chan_selected;
                ECOG=ft_selectdata(cfg,E);
            end
        case 'ecog'
            cfg=[];
            cfg.channel     ={'ecog*'};
            ECOG=ft_selectdata(cfg,E);
    end
    if cc==1
        cfg=[];
        cfg.label   = electrode.electrode;
        cfg.group   = electrode.connector;
        cfg.method  = 'CTAR';   % using trimmed average referencing
        cfg.percent = 50;       % percentage of 'extreme' channels in group to trim
        ECOG=bml_rereference(cfg,ECOG);
    end

    % dbs
    cfg=[];
    cfg.channel     ={'dbs_L*'};
    DBS=ft_selectdata(cfg,E);
    if strcmp(SUBJECT,"DBS3028")
        cfg=[];
        cfg.channel=['all', strcat('-', {'dbs_LS','dbs_L10'})];
        DBS=ft_selectdata(cfg,DBS);
    end
    DBS=DBS_referencing(DBS);

    % unify set of data
    cfg=[];
    E=ft_appenddata(cfg,ECOG,DBS);

    chan_selected=E.label;
    electrode=electrode(ismember(electrode.electrode,chan_selected),:);
    
    %% ft
    disp('*************************************************************')
    fprintf('Time Frequency Analysis for subject %s \n',SUBJECT)
    nTrials=numel(E.trial);
    criteria={'Volume','RT','Accuracy'};
    cond_definition={'positive','negative','difference'};
 
    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);

    nElec=numel(chan_selected);
    
    baselinetype='db'; % 'difference', 'percentage', 'db', 'zscore'
    for e=1:(numel(ECOG.label)+1)
        disp('---------------------------------------------------------')
        fprintf('%s %s',SUBJECT,chan_selected{e})
        
        % select single channel 
        cfg=[];
        cfg.channel     ={chan_selected{e}}; 
        Esingle=ft_selectdata(cfg,E);

        data=cell2mat(Esingle.trial');
        if sum(data(:)==0)/numel(data)>=0.7;disp('Electrode rejected, going to the next');continue;end

        % apply morlet wavelet
        cfg=[];
        cfg.foi     =TF_FOI;
        cfg.dt      =1/TF_RATE;
        cfg.toilim  =tf_epoch_t0;
        width=7;
        cfg.width   =width;
        Esingle_mor=bml_freqanalysis_power_wavelet(cfg,Esingle); 
        
        % xtime to plot (only original epoch)
        [~,aidx]=min(abs(Esingle_mor.time-a));
        [~,bidx]=min(abs(Esingle_mor.time-b));
        Esingle_mor.time=Esingle_mor.time(aidx:bidx);
        Esingle_mor.powspctrm=Esingle_mor.powspctrm(:,:,aidx:bidx);

        % Plot welch PSD
        f1=5;f2=40;
        freq_idx=Esingle_mor.freq>=f1 & Esingle_mor.freq<=f2;
        pow_toplot=squeeze(median(Esingle_mor.powspctrm,1,'omitnan')); 
        [~,aidx]=min(abs(Esingle_mor.time-x_syl1));[~,bidx]=min(abs(Esingle_mor.time-x_syl3_end));
        pow_speech=pow_toplot(freq_idx,aidx:bidx);
        [~,aidx]=min(abs(Esingle_mor.time-x_stim1));[~,bidx]=min(abs(Esingle_mor.time-x_stim3_end));
        pow_stim=pow_toplot(freq_idx,aidx:bidx);
        [~,aidx]=min(abs(Esingle_mor.time-mean(tf_baseline_t0(:,1),'omitnan')));[~,bidx]=min(abs(Esingle_mor.time-mean(tf_baseline_t0(:,2),'omitnan')));
        pow_base=pow_toplot(freq_idx,aidx:bidx);
        [~,aidx]=min(abs(Esingle_mor.time-tf_rebound_t0(1)));[~,bidx]=min(abs(Esingle_mor.time-tf_rebound_t0(2)));
        pow_rb=pow_toplot(freq_idx,aidx:bidx);
        toplot={pow_speech,pow_stim,pow_base,pow_rb};
        fig1=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
        subplot(4,3,1)
        for j=1:numel(toplot)
            w=mean(toplot{j},2,'omitnan');
            loglog(Esingle_mor.freq(freq_idx),w,'LineWidth',1);hold on
        end
        xline([beta(1) beta(1)],'Color','k','LineStyle','--','LineWidth',1);text(beta(1)+0.5,0.9*max(w),[num2str(beta(1))+" Hz"])
        xline([beta(2) beta(2)],'Color','k','LineStyle',':','LineWidth',0.5);text(beta(2)+0.5,0.9*max(w),[num2str(beta(2))+" Hz"])
        xline([beta(3) beta(3)],'Color','k','LineStyle','--','LineWidth',1);text(beta(3)+0.5,0.9*max(w),[num2str(beta(3))+" Hz"])
        xlabel('Frequency (Hz)')
        ylabel("PSD (mV/Hz)")
        legend('Speech','Stimulus','Baseline','Rebound')
        title("PSD for "+num2str(f1)+"-"+num2str(f2)+" Hz")         

        % normalize data
        cfg=[];
        cfg.baselinetype=baselinetype;
        cfg.baselinedef=tf_baseline_t0;
        Esingle_norm=baseline_normalization(cfg,Esingle_mor);
        
        % xtime to plot (only original epoch)
        [~,aidx]=min(abs(Esingle_norm.time-a));
        [~,bidx]=min(abs(Esingle_norm.time-b));
        Esingle_norm.time=Esingle_norm.time(aidx:bidx);
        Esingle_norm.powspctrm=Esingle_norm.powspctrm(:,aidx:bidx);

        % spectrogram plot after normalization
        subplot(4,3,[2,3])
        imagesc(Esingle_norm.time, Esingle_norm.freq, Esingle_norm.powspctrm);
        colormap(bluewhitered);colorbar %colormap(parula)
        xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',2)
        xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',2)
        xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',2)
        xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',2)
        yline([beta(1) beta(1)],'Color','c','LineWidth',1)
        yline([beta(3) beta(3)],'Color','c','LineWidth',1)
        xline([x_stim1 x_stim1],'Color','g','LineWidth',3)
        xline([x_stim2 x_stim2],'Color','g','LineWidth',1.5,'LineStyle',':')
        xline([x_stim3 x_stim3],'Color','g','LineWidth',1.5,'LineStyle',':')
        xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',2)
        xline([x_syl1 x_syl1],'Color','r','LineWidth',3)
        xline([x_syl2 x_syl2],'Color','r','LineWidth',1.5,'LineStyle',':')
        xline([x_syl3 x_syl3],'Color','r','LineWidth',1.5,'LineStyle',':')
        xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',2)
        ylim([4,250])
        set(gca,'YDir','normal') 
        set(gca,'YScale','log')
        [cmin,cmax]=caxis;
        title(['Wavelet [',baselinetype,' norm] spectrogram with cycles=',num2str(width), ' n trials=',num2str(nTrials)])


        % spectrogram per condition 
        for j=1:numel(criteria)
            switch criteria{j}
            case 'Volume'    
                data=prod_triplet.stim_volume(tokeepCod);
                positive=find(data); % high volume
                negative=find(~data); % low
            case 'RT'
                data=coding.syl1_onset(tokeepCod)-cue.stim1_ends(tokeepCue);
                positive=find(data>=median(data,'omitnan')); % high RT
                negative=find(data<median(data,'omitnan')); % low
            case 'Accuracy'
                try
                    data=prod_triplet.phonetic_ontarget(tokeepCod);
                    positive=find(data); % high accuracy
                    negative=find(~data); % low
                catch
                    continue
                end
            end
            for k=1:numel(cond_definition)
                if k<3
                    cfg=[];
                    if k==1;cfg.trials=positive';tit=strcat('N trials: ',num2str(length(positive)));
                    elseif k==2;cfg.trials=negative';tit=strcat('N trials: ',num2str(length(negative)));end
                    Econd=ft_selectdata(cfg,Esingle); 
                    
                    % tf analysis
                    cfg=[];
                    cfg.foi     =TF_FOI;
                    cfg.dt      =1/TF_RATE;
                    cfg.toilim  =tf_epoch_t0;
                    width=7;
                    cfg.width   =width;
                    Econd_mor=bml_freqanalysis_power_wavelet(cfg,Econd);
                    
                    % normalize data
                    cfg=[];
                    cfg.baselinetype=baselinetype;
                    cfg.baselinedef=tf_baseline_t0;
                    Econd_norm=baseline_normalization(cfg,Econd_mor);
                    
                    % xtime to plot (only original epoch)
                    [~,aidx]=min(abs(Econd_norm.time-a));
                    [~,bidx]=min(abs(Econd_norm.time-b));
                    Econd_norm.time=Econd_norm.time(aidx:bidx);
                    Econd_norm.powspctrm=Econd_norm.powspctrm(:,aidx:bidx);
                    
                    if k==1;Epos=Econd_norm;
                    elseif k==2;Eneg=Econd_norm;end
                    Etoplot=Econd_norm;
                else
                    Etoplot=Epos;Etoplot.powspctrm=Epos.powspctrm-Eneg.powspctrm;
                    tit='';
                end

                subplot(4,3,2*k+k+j)
                imagesc(Etoplot.time, Etoplot.freq, Etoplot.powspctrm);
                %colormap(bluewhitered);colorbar %colormap(parula)
                xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',2)
                xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',2)
                xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',2)
                xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',2)
                yline([beta(1) beta(1)],'Color','c','LineWidth',1)
                yline([beta(3) beta(3)],'Color','c','LineWidth',1)
                xline([x_stim1 x_stim1],'Color','g','LineWidth',3)
                xline([x_stim2 x_stim2],'Color','g','LineWidth',1.5,'LineStyle',':')
                xline([x_stim3 x_stim3],'Color','g','LineWidth',1.5,'LineStyle',':')
                xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',2)
                xline([x_syl1 x_syl1],'Color','r','LineWidth',3)
                xline([x_syl2 x_syl2],'Color','r','LineWidth',1.5,'LineStyle',':')
                xline([x_syl3 x_syl3],'Color','r','LineWidth',1.5,'LineStyle',':')
                xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',2)
                ylim([4,250])
                set(gca,'YDir','normal') 
                set(gca,'YScale','log')
                caxis([cmin, cmax])
                tit=[criteria{j},' ',cond_definition{k},' ',tit];
                title(tit)
            end
        end
        if e>numel(ECOG.label);tit=[strcat(SUBJECT,": ",strrep(string(Esingle_norm.label),'_',' ')),strcat('Locked in speech')];
        else
            try
                tit=[strcat(SUBJECT,": ",strrep(string(Esingle_norm.label),'_',' ')," ",string(electrode.HCPMMP1_area_fullname(e))),strcat('Locked in speech')];
            catch
                tit=[strcat(SUBJECT,": ",strrep(string(Esingle_norm.label),'_',' ')," [area not present]"),strcat('Locked in speech')];
            end
        end
        sgtitle(tit)
        [cmin,cmax]=caxis;
        if cc==1;name=[strcat('images/conditions CTAR/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_normdb')];
        elseif cc==2;name=[strcat('images/conditions/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_normdb')];
        end
        saveas(fig1,strcat(name,'.png'))
        saveas(fig1,strcat(name,'.fig'))
        
        close all
    end

end
end
%  pcolor(Esingle_norm.time, Esingle_norm.freq,pow)
%  set(gca,'yscale','log')      




