close all;clear all;clc;
%% 
%%% this script get raw data and 
% 1) select session id with DBS LFP data
% 2) preprocess 
% 3) epoch [syl1, stim1]
% 4) filter out artifact from annot  
% 5) save

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
for i=1:numel(SUBJECTS)
    if i==6;continue;end
    %open a subject dir
    tic
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
    
    % define epoch and baseline for ft
    tf_epoch_t0=E.time{:};tf_epoch_t0=[tf_epoch_t0(1) tf_epoch_t0(end)];
    % tf_epoch_t0=[info_epoch.starts-center info_epoch.ends-center];
    % tf_epoch_t0=[info_epoch.starts(1) info_epoch.ends(1)];
    tf_epoch_t0=round(tf_epoch_t0 .* TF_RATE) ./ TF_RATE;

    x=center - cue.stim1_starts(tokeepCue); %TRIAL SPECIFIC
    % x=median(coding.syl1_onset - cue.stim1_starts,'omitnan'); %MEDIATING
    tf_baseline_t0=[-(x+0.9) -(x+0.1)];
    tf_baseline_t0=round(tf_baseline_t0 .* TF_RATE) ./ TF_RATE;

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
        % save(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_rebound_clean_ref_globt.mat"), 'E_ctar')
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
    electrode=electrode(ismember(electrode.electrode,ECOG.label),:);
            
    %% ft
    disp('*************************************************************')
    fprintf('Time Frequency Analysis for subject %s \n',SUBJECT)
    nTrials=numel(E.trial);
    powtab=table();
    powtabp=table();

    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);

    nElec=numel(chan_selected);
    
    for e=1:nElec
        disp('---------------------------------------------------------')
        fprintf('%s %s',SUBJECT,chan_selected{e})
        
        % select single channel and continue if it's a rejected channel
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
        
        % normalize data db
        cfg=[];
        cfg.baselinetype='db';
        cfg.baselinedef=tf_baseline_t0;
        Esingle_norm=baseline_normalization(cfg,Esingle_mor);

        % xtime to plot (only original epoch)
        [~,aidx]=min(abs(Esingle_norm.time-a));
        [~,bidx]=min(abs(Esingle_norm.time-b));
        Esingle_norm.time=Esingle_norm.time(aidx:bidx);
        Esingle_norm.powspctrm=Esingle_norm.powspctrm(:,aidx:bidx);

     
        % normalize data percentage
        cfg=[];
        cfg.baselinetype='percentage';
        cfg.baselinedef=tf_baseline_t0;
        Esingle_normp=baseline_normalization(cfg,Esingle_mor);

        % xtime to plot (only original epoch)
        [~,aidx]=min(abs(Esingle_normp.time-a));
        [~,bidx]=min(abs(Esingle_normp.time-b));
        Esingle_normp.time=Esingle_normp.time(aidx:bidx);
        Esingle_normp.powspctrm=Esingle_normp.powspctrm(:,aidx:bidx);

    
        % mean power in beta during task and rebound
        pow_beta=Esingle_norm.powspctrm( Esingle_norm.freq>=beta(1)&Esingle_norm.freq<=beta(3) ,:);
        pow_highbeta=Esingle_norm.powspctrm( Esingle_norm.freq>=beta(2)&Esingle_norm.freq<=beta(3) ,:);
        pow_lowbeta=Esingle_norm.powspctrm( Esingle_norm.freq>=beta(1)&Esingle_norm.freq<=beta(2) ,:);
        %powtab.id(e)=electrode.id(strcmp(electrode.electrode,chan_selected{e}));
        powtab.id(e)=e;
        powtab.label(e)=chan_selected(e);
            %baseline 
        [~,idx1]=min(abs(Esingle_norm.time-mean(tf_baseline_t0(:,1),'omitnan')));
        [~,idx2]=min(abs(Esingle_norm.time-mean(tf_baseline_t0(:,2),'omitnan')));
        powtab.base_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtab.base_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtab.base_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % stim
        [~,idx1]=min(abs(Esingle_norm.time-x_stim1));
        [~,idx2]=min(abs(Esingle_norm.time-x_stim3_end));
        powtab.stim_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtab.stim_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtab.stim_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % prespeech
        [~,idx1]=min(abs(Esingle_norm.time-x_stim3_end));
        [~,idx2]=min(abs(Esingle_norm.time-x_syl1));
        powtab.prespeech_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtab.prespeech_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtab.prespeech_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % speech
        [~,idx1]=min(abs(Esingle_norm.time-x_syl1));
        [~,idx2]=min(abs(Esingle_norm.time-x_syl3_end));
        powtab.speech_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtab.speech_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtab.speech_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');

            % rebound
        [~,idx1]=min(abs(Esingle_norm.time-tf_rebound_t0(1)));
        [~,idx2]=min(abs(Esingle_norm.time-tf_rebound_t0(2)));
        powtab.rebound_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtab.rebound_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtab.rebound_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');  


        % mean power in beta during task and rebound PERCENTAGE
        pow_beta=Esingle_normp.powspctrm( Esingle_normp.freq>=beta(1)&Esingle_normp.freq<=beta(3) ,:);
        pow_highbeta=Esingle_normp.powspctrm( Esingle_normp.freq>=beta(2)&Esingle_normp.freq<=beta(3) ,:);
        pow_lowbeta=Esingle_normp.powspctrm( Esingle_normp.freq>=beta(1)&Esingle_normp.freq<=beta(2) ,:);
        powtabp.id(e)=e;
        powtabp.label(e)=chan_selected(e);
            %baseline 
        [~,idx1]=min(abs(Esingle_normp.time-mean(tf_baseline_t0(:,1),'omitnan')));
        [~,idx2]=min(abs(Esingle_normp.time-mean(tf_baseline_t0(:,2),'omitnan')));
        powtabp.base_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtabp.base_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtabp.base_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % stim
        [~,idx1]=min(abs(Esingle_normp.time-x_stim1));
        [~,idx2]=min(abs(Esingle_normp.time-x_stim3_end));
        powtabp.stim_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtabp.stim_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtabp.stim_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % prespeech
        [~,idx1]=min(abs(Esingle_normp.time-x_stim3_end));
        [~,idx2]=min(abs(Esingle_normp.time-x_syl1));
        powtabp.prespeech_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtabp.prespeech_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtabp.prespeech_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % speech
        [~,idx1]=min(abs(Esingle_normp.time-x_syl1));
        [~,idx2]=min(abs(Esingle_normp.time-x_syl3_end));
        powtabp.speech_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtabp.speech_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtabp.speech_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan');
            % rebound post speech
        [~,idx1]=min(abs(Esingle_normp.time-tf_rebound_t0(1)));
        [~,idx2]=min(abs(Esingle_normp.time-tf_rebound_t0(2)));
        powtabp.rebound_beta(e)=mean(pow_beta(:, idx1:idx2),'all','omitnan');
        powtabp.rebound_highbeta(e)=mean(pow_highbeta(:, idx1:idx2),'all','omitnan');
        powtabp.rebound_lowbeta(e)=mean(pow_lowbeta(:, idx1:idx2),'all','omitnan'); 
    end
    powtab=powtab(~powtab.id==0, :);
    powtabp=powtabp(~powtabp.id==0, :);
    powtab.id=(1:height(powtab))';
    powtabp.id=(1:height(powtabp))';

    if cc==1
        filename = strcat(SUBJECT,'_speech_pow_db.txt');
        writetable(powtab, strcat('annot/general CTAR/',filename), 'Delimiter', '\t');
        filename = strcat(SUBJECT,'_speech_pow_percentage.txt');
        writetable(powtabp, strcat('annot/general CTAR/',filename), 'Delimiter', '\t');
    elseif cc==2
        filename = strcat(SUBJECT,'_speech_pow_db.txt');
        writetable(powtab, strcat('annot/general/',filename), 'Delimiter', '\t');
        filename = strcat(SUBJECT,'_speech_pow_percentage.txt');
        writetable(powtabp, strcat('annot/general/',filename), 'Delimiter', '\t');
    end
    disp(toc)
end
end
disp('The end :)')
      





