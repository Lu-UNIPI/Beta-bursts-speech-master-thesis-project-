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

% freq analysis settings
TF_RATE = 20; %Hz Sampling rate of time frequency plot
TF_FOI = round(10.^(0.30:0.05:2.4),2,'signif');
TF_BASELINE_WIDTH = 0.5;
beta=[11,22,36];
percTh=75;
deltab=1.2;

% FOOOF settings
settings_fooof = struct();  % Use defaults
f_range_fooof = [5, 50];
f_range_bycycle = [10, 40];

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for cc=1
for i=5%1:numel(SUBJECTS) 
    %open a subject dir
    if i==6;continue;end
    %if i==14;continue;end
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
        case {'3012','3014'};tokeepCod=tokeepCue;
        case {'3024'};tokeepCod=tokeepCod & tokeepCue;tokeepCue=tokeepCod;
        case {'3021'};idx=find(~tokeepCue)+1;idx=idx(1:2);tokeepCod(idx)=0;
    end

    % load data
    load(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_syl.mat"));a=-(A+LEFT);b=B+RIGHT;center=coding.syl1_onset(tokeepCod);
    %load(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_stim.mat"));a=-(IT+LEFT);b=A+RIGHT;center=cue.stim1_starts(tokeepCue);
    fsample=E.fsample;

    % define epoch, baseline and rebound
    tf_epoch_t0=E.time{:};tf_epoch_t0=[tf_epoch_t0(1) tf_epoch_t0(end)];
    tf_epoch_t0=round(tf_epoch_t0 .* TF_RATE) ./ TF_RATE;

    x=center - cue.stim1_starts(tokeepCue); %TRIAL SPECIFIC
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
    rt=1;
    
    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);

    nElec=numel(chan_selected); 
    ifreqs=nan(1,nElec);
    for e=31%1:nElec
        disp('---------------------------------------------------------')
        fprintf('%s %s',SUBJECT,chan_selected{e})
        pow_beta={};mask_beta={};
        pow_fooof={};mask_fooof={};
        pow_time={};mask_time={};
        deltaf=nan(nTrials,1);
        % select single channel and continue if it's a rejected channel
        cfg=[];
        cfg.channel     ={chan_selected{e}}; 
        Esingle0=ft_selectdata(cfg,E);
        
        data0=cell2mat(Esingle0.trial');
        if sum(data0(:)==0)/numel(data0)>=0.7;disp('Electrode rejected, going to the next');continue;end
        [~,aidx]=min(abs(E.time{1}-a));
        [~,bidx]=min(abs(E.time{1}-b));
        Etime=E.time{1}(aidx:bidx);
        data0=data0(:,aidx:bidx);
        deltat=prctile(abs(data0),90,'all');   

        window_pxx=256;noverlap_pxx=128;nfft_pxx=512;fs_pxx=Esingle0.fsample;
        [pxx,freq_pxx]=cellfun(@(x) pwelch(x,window_pxx,noverlap_pxx,nfft_pxx,fs_pxx),Esingle0.trial,'UniformOutput',false);
        freq_pxx=freq_pxx{1};freq_pxx=freq_pxx';
        pxx=cell2mat(pxx);pxx=mean(pxx,2,'omitnan');

        % iFreq from fooof
        try
            fooof_results=fooof(freq_pxx, pxx, f_range_fooof, settings_fooof);
            % NB data must not contain values <=0, or log(data) will result in NaN or Inf
            if ~isempty(fooof_results.peak_params) % case there are results in fooof
                fooof_peak_params=fooof_results.peak_params( fooof_results.peak_params(:,1)>=f_range_bycycle(1) & fooof_results.peak_params(:,1)<=f_range_bycycle(2) , :);
                [~,max_peak]=max(fooof_peak_params(:,2).*fooof_peak_params(:,3));
                iFreq=fooof_peak_params(max_peak,1);
                if isempty(iFreq) % case fooof results are outside freqlimits
                    mask=freq_pxx<f_range_bycycle(2) & freq_pxx>f_range_bycycle(1);
                    [~,max_peak]=max(abs(pxx(mask)));
                    id=find(mask);
                    iFreq=freq_pxx(id(1)-1+max_peak);
                end
            else
                mask=freq_pxx<f_range_bycycle(2) & freq_pxx>f_range_bycycle(1);
                [~,max_peak]=max(abs(pxx(mask)));
                id=find(mask);
                iFreq=freq_pxx(id(1)-1+max_peak);
            end
            ifreqs(e)=iFreq;
        catch
            mask=freq_pxx<f_range_bycycle(2) & freq_pxx>f_range_bycycle(1);
            [~,max_peak]=max(abs(pxx(mask)));
            id=find(mask);
            iFreq=freq_pxx(id(1)-1+max_peak);
            ifreqs(e)=iFreq;
        end

        % get data
        for j=1:nTrials % for every trial
            cfg=[];
            cfg.trials=[1:numel(Esingle0.trial)]==j;
            Esingle=ft_selectdata(cfg,Esingle0);
            data=cell2mat(Esingle.trial');
            if sum(data(:)==0)/numel(data)>=0.95
                disp('Trial rejected, going to the next');
                if j==nTrials;pow_beta{j}=[];pow_fooof{j}=[];ifreqs(j)=NaN;end
                continue
            end

            % apply morlet wavelet
            cfg=[];
            cfg.foi     =TF_FOI;
            cfg.dt      =1/TF_RATE;
            cfg.toilim  =tf_epoch_t0;
            width=7;
            cfg.width   =width;
            Esingle_mor=bml_freqanalysis_power_wavelet(cfg,Esingle); 
            
            
            % normalize data zscore
            cfg=[];
            cfg.baselinetype='zscore';
            Esingle_norm=baseline_normalization(cfg,Esingle_mor);
            % xtime to plot (only original epoch)
            [~,aidx]=min(abs(Esingle_norm.time-a));
            [~,bidx]=min(abs(Esingle_norm.time-b));
            Esingle_norm.time=Esingle_norm.time(aidx:bidx);
            Esingle_norm.powspctrm=Esingle_norm.powspctrm(:,aidx:bidx);

            % 1) get bursts by freq domain general
            mask=Esingle_norm.freq<35 & Esingle_norm.freq>12;
            powBeta=mean(Esingle_norm.powspctrm(mask,:),'omitnan');
            pow_beta{j}=powBeta;

            % 2) get bursts by freq domain iFreq from fooof
            % iFreq from fooof
            % NB data must not contain values <=0, or log(data) will result in NaN or Inf
            maskf=Esingle_norm.freq<ifreqs(e)+2.5 & Esingle_norm.freq>ifreqs(e)-2.5;
            if sum(maskf)==1;powFooof=Esingle_norm.powspctrm(maskf,:);
            else;powFooof=mean(Esingle_norm.powspctrm(maskf,:),'omitnan');
            end            
            pow_fooof{j}=powFooof;
            deltaf(j)=prctile(mean(Esingle_norm.powspctrm,2,'omitnan'),90,'all');  
        end
        deltaf=mean(deltaf,'omitnan');
        
        %% test frequency parameters
        percentiles=20:5:95;
        n_beta=nan(length(percentiles),nTrials);n_fooof=nan(length(percentiles),nTrials);
        d_beta=nan(length(percentiles),nTrials);d_fooof=nan(length(percentiles),nTrials);
        a_beta=nan(length(percentiles),nTrials);a_fooof=nan(length(percentiles),nTrials);
        mp_beta=nan(length(percentiles),nTrials);mp_fooof=nan(length(percentiles),nTrials);
        mask_beta={};mask_fooof={};
        beta_baseline=cell(1,4);beta_stimulus=cell(1,4);beta_speech=cell(1,4);beta_rebound=cell(1,4);
        fooof_baseline=cell(1,4);fooof_stimulus=cell(1,4);fooof_speech=cell(1,4);fooof_rebound=cell(1,4);

        [~,bl_idx]=min(abs(Esingle_norm.time-mean(tf_baseline_t0(:,1),'omitnan')));
        [~,br_idx]=min(abs(Esingle_norm.time-mean(tf_baseline_t0(:,2),'omitnan')));
        [~,stl_idx]=min(abs(Esingle_norm.time-x_stim1));
        [~,str_idx]=min(abs(Esingle_norm.time-x_stim3_end));
        [~,syl_idx]=min(abs(Esingle_norm.time-x_syl1));
        [~,syr_idx]=min(abs(Esingle_norm.time-x_syl3_end));
        [~,rl_idx]=min(abs(Esingle_norm.time-tf_rebound_t0(1)));
        [~,rr_idx]=min(abs(Esingle_norm.time-tf_rebound_t0(2)));

        for p=percentiles
            percBeta=prctile(cell2mat(pow_beta'),p,'all');
            percFooof=prctile(cell2mat(pow_fooof'),p,'all');
            for j=1:nTrials
                if isempty(pow_beta{j});continue;end

                % 1) get bursts by freq domain general
                cfg=[];
                cfg.threshold=percBeta;
                cfg.fs=TF_RATE;
                cfg.domain='frequency';
                maskTh=bursts_consistence(cfg,pow_beta{j});
                mask_beta{j}=maskTh;
                pieces=bwconncomp(maskTh);
                n_beta(find(percentiles==p),j)=pieces.NumObjects;
                duration=nan(1,pieces.NumObjects);amplitude=nan(1,pieces.NumObjects);maxpow=nan(1,pieces.NumObjects);
                for ii=1:pieces.NumObjects
                    idx=pieces.PixelIdxList{ii};
                    duration(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));
                    amplitude(ii)=1/TF_RATE*trapz(pow_beta{j}(idx))/duration(ii);
                    maxpow(ii)=max(pow_beta{j}(idx));
                end
                d_beta(find(percentiles==p),j)=mean(duration,'omitnan');
                a_beta(find(percentiles==p),j)=mean(amplitude,'omitnan');
                mp_beta(find(percentiles==p),j)=mean(maxpow,'omitnan');
                
                cfg=[];
                cfg.limits=[bl_idx,br_idx,stl_idx,str_idx,syl_idx,syr_idx,rl_idx,rr_idx];
                pieces=window_belonging(cfg,pieces);
                beta_baseline{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'baseline'));
                beta_stimulus{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'stimulus'));
                beta_speech{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'speech'));
                beta_rebound{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'rebound'));
                duration_b=nan(1,pieces.NumObjects);duration_st=nan(1,pieces.NumObjects);duration_sy=nan(1,pieces.NumObjects);duration_r=nan(1,pieces.NumObjects);
                amplitude_b=nan(1,pieces.NumObjects);amplitude_st=nan(1,pieces.NumObjects);amplitude_sy=nan(1,pieces.NumObjects);amplitude_r=nan(1,pieces.NumObjects);
                maxpow_b=nan(1,pieces.NumObjects);maxpow_st=nan(1,pieces.NumObjects);maxpow_sy=nan(1,pieces.NumObjects);maxpow_r=nan(1,pieces.NumObjects);
                for ii=1:pieces.NumObjects
                    switch pieces.area{ii}
                        case 'baseline'
                            idx=pieces.PixelIdxList{ii};
                            duration_b(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));
                            amplitude_b(ii)=1/TF_RATE*trapz(pow_beta{j}(idx))/duration_b(ii);
                            maxpow_b(ii)=max(pow_beta{j}(idx));
                        case 'stimulus'
                            idx=pieces.PixelIdxList{ii};
                            duration_st(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1)); 
                            amplitude_st(ii)=1/TF_RATE*trapz(pow_beta{j}(idx))/duration_st(ii);
                            maxpow_st(ii)=max(pow_beta{j}(idx));
                        case 'speech'
                            idx=pieces.PixelIdxList{ii};
                            duration_sy(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));    
                            amplitude_sy(ii)=1/TF_RATE*trapz(pow_beta{j}(idx))/duration_sy(ii);
                            maxpow_sy(ii)=max(pow_beta{j}(idx));
                        case 'rebound'
                            idx=pieces.PixelIdxList{ii};
                            duration_r(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));  
                            amplitude_r(ii)=1/TF_RATE*trapz(pow_beta{j}(idx))/duration_r(ii);
                            maxpow_r(ii)=max(pow_beta{j}(idx));
                    end
                end
                beta_baseline{2}(find(percentiles==p),j)=mean(duration_b,'omitnan');
                beta_stimulus{2}(find(percentiles==p),j)=mean(duration_st,'omitnan');
                beta_speech{2}(find(percentiles==p),j)=mean(duration_sy,'omitnan');
                beta_rebound{2}(find(percentiles==p),j)=mean(duration_r,'omitnan');
                beta_baseline{3}(find(percentiles==p),j)=mean(amplitude_b,'omitnan');
                beta_stimulus{3}(find(percentiles==p),j)=mean(amplitude_st,'omitnan');
                beta_speech{3}(find(percentiles==p),j)=mean(amplitude_sy,'omitnan');
                beta_rebound{3}(find(percentiles==p),j)=mean(amplitude_r,'omitnan');      
                beta_baseline{4}(find(percentiles==p),j)=mean(maxpow_b,'omitnan');
                beta_stimulus{4}(find(percentiles==p),j)=mean(maxpow_st,'omitnan');
                beta_speech{4}(find(percentiles==p),j)=mean(maxpow_sy,'omitnan');
                beta_rebound{4}(find(percentiles==p),j)=mean(maxpow_r,'omitnan');  

                % 2) get bursts by freq domain iFreq from fooof
                cfg=[];
                cfg.threshold=percFooof;
                cfg.fs=TF_RATE;
                cfg.domain='frequency';
                maskTh=bursts_consistence(cfg,pow_fooof{j});
                mask_fooof{j}=maskTh;
                pieces=bwconncomp(maskTh);
                n_fooof(find(percentiles==p),j)=pieces.NumObjects;
                duration=nan(1,pieces.NumObjects);amplitude=nan(1,pieces.NumObjects);maxpow=nan(1,pieces.NumObjects);
                for ii=1:pieces.NumObjects
                    idx=pieces.PixelIdxList{ii};
                    duration(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));
                    amplitude(ii)=1/TF_RATE*trapz(pow_fooof{j}(idx))/duration(ii);
                    maxpow(ii)=max(pow_fooof{j}(idx));
                end
                d_fooof(find(percentiles==p),j)=mean(duration,'omitnan');
                a_fooof(find(percentiles==p),j)=mean(amplitude,'omitnan');
                mp_fooof(find(percentiles==p),j)=mean(maxpow,'omitnan');

                cfg=[];
                cfg.limits=[bl_idx,br_idx,stl_idx,str_idx,syl_idx,syr_idx,rl_idx,rr_idx];
                pieces=window_belonging(cfg,pieces);
                fooof_baseline{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'baseline'));
                fooof_stimulus{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'stimulus'));
                fooof_speech{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'speech'));
                fooof_rebound{1}(find(percentiles==p),j)=sum(strcmp(pieces.area,'rebound'));
                duration_b=nan(1,pieces.NumObjects);duration_st=nan(1,pieces.NumObjects);duration_sy=nan(1,pieces.NumObjects);duration_r=nan(1,pieces.NumObjects);
                amplitude_b=nan(1,pieces.NumObjects);amplitude_st=nan(1,pieces.NumObjects);amplitude_sy=nan(1,pieces.NumObjects);amplitude_r=nan(1,pieces.NumObjects);
                maxpow_b=nan(1,pieces.NumObjects);maxpow_st=nan(1,pieces.NumObjects);maxpow_sy=nan(1,pieces.NumObjects);maxpow_r=nan(1,pieces.NumObjects);
                for ii=1:pieces.NumObjects
                    switch pieces.area{ii}
                        case 'baseline'
                            idx=pieces.PixelIdxList{ii};
                            duration_b(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));
                            amplitude_b(ii)=1/TF_RATE*trapz(pow_fooof{j}(idx))/duration_b(ii);
                            maxpow_b(ii)=max(pow_fooof{j}(idx));
                        case 'stimulus'
                            idx=pieces.PixelIdxList{ii};
                            duration_st(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1)); 
                            amplitude_st(ii)=1/TF_RATE*trapz(pow_fooof{j}(idx))/duration_st(ii);
                            maxpow_st(ii)=max(pow_fooof{j}(idx));
                        case 'speech'
                            idx=pieces.PixelIdxList{ii};
                            duration_sy(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));    
                            amplitude_sy(ii)=1/TF_RATE*trapz(pow_fooof{j}(idx))/duration_sy(ii);
                            maxpow_sy(ii)=max(pow_fooof{j}(idx));
                        case 'rebound'
                            idx=pieces.PixelIdxList{ii};
                            duration_r(ii)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));  
                            amplitude_r(ii)=1/TF_RATE*trapz(pow_fooof{j}(idx))/duration_r(ii);
                            maxpow_r(ii)=max(pow_fooof{j}(idx));
                    end
                end
                fooof_baseline{2}(find(percentiles==p),j)=mean(duration_b,'omitnan');
                fooof_stimulus{2}(find(percentiles==p),j)=mean(duration_st,'omitnan');
                fooof_speech{2}(find(percentiles==p),j)=mean(duration_sy,'omitnan');
                fooof_rebound{2}(find(percentiles==p),j)=mean(duration_r,'omitnan');
                fooof_baseline{3}(find(percentiles==p),j)=mean(amplitude_b,'omitnan');
                fooof_stimulus{3}(find(percentiles==p),j)=mean(amplitude_st,'omitnan');
                fooof_speech{3}(find(percentiles==p),j)=mean(amplitude_sy,'omitnan');
                fooof_rebound{3}(find(percentiles==p),j)=mean(amplitude_r,'omitnan');
                fooof_baseline{4}(find(percentiles==p),j)=mean(maxpow_b,'omitnan');
                fooof_stimulus{4}(find(percentiles==p),j)=mean(maxpow_st,'omitnan');
                fooof_speech{4}(find(percentiles==p),j)=mean(maxpow_sy,'omitnan');
                fooof_rebound{4}(find(percentiles==p),j)=mean(maxpow_r,'omitnan');
            end
            mask_beta=mask_beta(~cellfun(@isempty, mask_beta));
            m_beta{find(percentiles==p)}=mean(cell2mat(mask_beta'),'omitnan');

            mask_fooof=mask_fooof(~cellfun(@isempty, mask_fooof));
            m_fooof{find(percentiles==p)}=mean(cell2mat(mask_fooof'),'omitnan');
        end
       
        
        %% FIGURE: parameters general
        fig1=figure('units','normalized','outerposition',[0.01 0.01 0.99 0.99]);
        subplot(3,8,1);
        sig=mean(n_beta,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('mean n bursts')
        title('Mean n bursts')
        subplot(3,8,2)
        sig=mean(d_beta,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('Burst duration')
        title({'General beta', '','Burst duration'})
        subplot(3,8,3)
        sig=mean(a_beta,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('Burst amplitude')
        title('Burst ampitude')
        subplot(3,8,4)
        sig=mean(mp_beta,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('Burst max pow')
        title('Burst max pow')
        subplot(3,8,[9,10,11,12])
        imagesc([Esingle_norm.time(1) Esingle_norm.time(end)],[percentiles(1) percentiles(end)],cell2mat(m_beta'));colorbar
        xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1)
        xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1)
        xline([x_stim1 x_stim1],'Color','g','LineWidth',1)
        xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1)
        xline([x_syl1 x_syl1],'Color','r','LineWidth',1)
        xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1)
        yline(75,'Color','y','LineWidth',2)
        xlabel('time [s]')
        ylabel('Percentile')
        title('General beta - burst probability')

        h1=subplot(3,8,[17,18,19,20]);
        pos = get(h1, 'Position');
        new_pos = [pos(1) pos(2) 0.9*pos(3) 0.92*pos(4)]; 
        set(h1, 'Position', new_pos);
        axis tight
        plot(Esingle_norm.time,m_beta{find(percentiles==75)});
        ylabel('Bursts prob. 75th')
        yyaxis right
        plot(Esingle_norm.time,median(cell2mat(pow_beta')));hold on
        xlim([mean(tf_baseline_t0(:,1),'omitnan')-0.5 tf_rebound_t0(2)+0.5])
        xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1)
        xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1)
        xline([x_stim1 x_stim1],'Color','g','LineWidth',1)
        xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1)
        xline([x_syl1 x_syl1],'Color','r','LineWidth',1)
        xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1)
        legend('Burst prob at 75 perc','Beta power')
        xlabel('time [s]')
        ylabel('Beta power')
        title('General beta - 75 percentile')

        subplot(3,8,5);
        sig=mean(n_fooof,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('mean n bursts')
        title('Mean n bursts')
        subplot(3,8,6)
        sig=mean(d_fooof,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('Burst duration')
        title({'Fooof iFreq', '','Burst duration'})
        subplot(3,8,7)
        sig=mean(a_fooof,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('Burst amplitude')
        title('Burst ampitude')
        subplot(3,8,8)
        sig=mean(mp_fooof,2,'omitnan');plot(percentiles,sig)
        xline(75)
        ylim([0,max(sig)*1.1])
        xlabel('Percentile')
        ylabel('Burst max pow')
        title('Burst max pow')
        subplot(3,8,[13,14,15,16])
        imagesc([Esingle_norm.time(1) Esingle_norm.time(end)],[percentiles(1) percentiles(end)],cell2mat(m_fooof'));colormap(linspecer);colorbar
        xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1)
        xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1)
        xline([x_stim1 x_stim1],'Color','g','LineWidth',1)
        xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1)
        xline([x_syl1 x_syl1],'Color','r','LineWidth',1)
        xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1)
        yline(75,'Color','y','LineWidth',2)
        xlabel('time [s]')
        ylabel('Percentile')
        title('Fooof iFreq - burst probability')

        h2=subplot(3,8,[21,22,23,24]);
        pos = get(h2, 'Position');
        new_pos = [pos(1) pos(2) 0.9*pos(3) 0.92*pos(4)]; 
        set(h2, 'Position', new_pos);
        plot(Esingle_norm.time,m_fooof{find(percentiles==75)});
        ylabel('Bursts prob. 75th')
        yyaxis right
        plot(Esingle_norm.time,median(cell2mat(pow_fooof'),'omitnan'));hold on
        xlim([mean(tf_baseline_t0(:,1),'omitnan')-0.5 tf_rebound_t0(2)+0.5])
        xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1)
        xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1)
        xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1)
        xline([x_stim1 x_stim1],'Color','g','LineWidth',1)
        xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
        xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1)
        xline([x_syl1 x_syl1],'Color','r','LineWidth',1)
        xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
        xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1)
        legend('Burst prob at 75 perc','Fooof power')
        xlabel('time [s]')
        ylabel('iFreq power')
        %ylabel('Power')
        title('Fooof iFreq - 75 percentile')
        sgtitle([SUBJECT+": "+strrep(chan_selected{e},'_','')+' '+'frequency parameter comparison per percentile'])
        name=[strcat('images/bursts CTAR/frequency parameters//',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_freq_params_percentile')];
        saveas(fig1,strcat(name,'.png'))
        saveas(fig1,strcat(name,'.fig'))

        %% FIGURE: Parameters in windows
        fig2=figure('units','normalized','outerposition',[0.01 0.01 0.99 0.99]);
        n_plot=[1,3,5,7];
        for jj=1:length(n_plot)
            maxx=max([max(mean(beta_baseline{jj},2,'omitnan'),[],'all'),max(mean(beta_stimulus{jj},2,'omitnan'),[],'all'),max(mean(beta_speech{jj},2,'omitnan'),[],'all'),max(mean(beta_rebound{jj},2,'omitnan'),[],'all'),...
                max(mean(fooof_baseline{jj},2,'omitnan'),[],'all'),max(mean(fooof_stimulus{jj},2,'omitnan'),[],'all'),max(mean(fooof_speech{jj},2,'omitnan'),[],'all'),max(mean(fooof_rebound{jj},2,'omitnan'),[],'all')]);
            if jj==1;max_n=maxx;
            elseif jj==2;max_d=maxx;
            elseif jj==3; max_a=maxx;
            elseif jj==4;max_p=maxx;
            end
        end
        for jj=1:length(n_plot)
            subplot(4,2,n_plot(jj))
            plot(percentiles,mean(beta_baseline{jj},2,'omitnan'));hold on
            plot(percentiles,mean(beta_stimulus{jj},2,'omitnan'));
            plot(percentiles,mean(beta_speech{jj},2,'omitnan'));
            plot(percentiles,mean(beta_rebound{jj},2,'omitnan'));
            xline(75)
            legend('baseline','stimulus','speech','rebound','Location','northwest')
            xlabel('Percentile')
            if jj==1;ylim([0,max_n]);ylabel('Mean n bursts');title({'General beta','','Mean n bursts'})
            elseif jj==2;ylim([0,max_d]);ylabel('Burst duration');title('Burst duration')
            elseif jj==3;ylim([0,max_a]);ylabel('Burst norm amplitude');title('Burst norm amplitude')
            elseif jj==4;ylim([0,max_p]);ylabel('Burst max pow');title('Burst max pow')
            end
        end
        n_plot=[2,4,6,8];
        for jj=1:length(n_plot)
            subplot(4,2,n_plot(jj))
            plot(percentiles,mean(fooof_baseline{jj},2,'omitnan'));hold on
            plot(percentiles,mean(fooof_stimulus{jj},2,'omitnan'));
            plot(percentiles,mean(fooof_speech{jj},2,'omitnan'));
            plot(percentiles,mean(fooof_rebound{jj},2,'omitnan'));
            xline(75)
            legend('baseline','stimulus','speech','rebound','Location','northwest')
            xlabel('Percentile')
            if jj==1;ylim([0,max_n]);ylabel('Mean n bursts');title({'Fooof iFreq','','Mean n bursts'})
            elseif jj==2;ylim([0,max_d]);ylabel('Burst duration');title('Burst duration')
            elseif jj==3;ylim([0,max_a]);ylabel('Burst norm amplitude');title('Burst norm amplitude')
            elseif jj==4;ylim([0,max_p]);ylabel('Burst max pow');title('Burst max pow')
            end
        end
        sgtitle([SUBJECT+": "+strrep(chan_selected{e},'_',' ')+' '+'frequency parameter comparison per window'])
        name=[strcat('images/bursts CTAR/frequency parameters/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_freq_params_in_time')];        
        saveas(fig2,strcat(name,'.png'))
        saveas(fig2,strcat(name,'.fig'))
        close all
    end
end
end