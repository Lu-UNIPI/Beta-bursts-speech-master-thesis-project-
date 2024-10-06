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

% BYCYCLE settings
center_extrema = 'trough';
switch center_extrema
    case 'trough';first_sample='sample_last_peak';last_sample='sample_next_peak';
    case 'peak';first_sample='sample_last_trough';last_sample='sample_next_trough';
end
burst_method = 'cycles';
settings_byclycle = struct();  
settings_bycycle.amp_fraction=0.4;
settings_bycycle.amp_consistency=0.35;
settings_bycycle.period_consistency=0.55;
settings_bycycle.monotonicity=0.6;
f_range_bycycle = [10, 40];
return_model = false;

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for cc=1
for i=15:16%:numel(SUBJECTS)
    %open a subject dir
    if i==6;continue;end
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
    bb=1;tt=1;ff=1;
    annot_beta=table();
    annot_fooof=table();
    annot_cycle=table();
    
    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);
%%
    nElec=numel(chan_selected); 
    ifreqs=nan(1,nElec);
    
    for e=1:nElec
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
            disp(strcat('Getting frequency data for trial n',string(j)))
            cfg=[];
            cfg.trials=[1:numel(Esingle0.trial)]==j;
            Esingle=ft_selectdata(cfg,Esingle0);
            data=cell2mat(Esingle.trial');
            if sum(data(:)==0)/numel(data)>=0.7
                pow_beta{j}=[];pow_fooof{j}=[];
                %if j==nTrials;pow_beta{j}=[];pow_fooof{j}=[];end
                disp('Trial rejected, going to the next');
                continue;
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
            mask=Esingle_norm.freq<beta(3) & Esingle_norm.freq>beta(1);
            powBeta=mean(Esingle_norm.powspctrm(mask,:),'omitnan');
            pow_beta{j}=powBeta;

            % 2) get bursts by freq domain iFreq from fooof
            maskf=Esingle_norm.freq<ifreqs(e)+2.5 & Esingle_norm.freq>ifreqs(e)-2.5;
            if sum(maskf)==1;powFooof=Esingle_norm.powspctrm(maskf,:);
            else;powFooof=mean(Esingle_norm.powspctrm(maskf,:),'omitnan');
            end            
            pow_fooof{j}=powFooof;
            deltaf(j)=prctile(mean(Esingle_norm.powspctrm,2,'omitnan'),90,'all');    
        end
        deltaf=mean(deltaf,'omitnan');
        percBeta=prctile(cell2mat(pow_beta'),75,'all');
        percFooof=prctile(cell2mat(pow_fooof'),75,'all');

        %% upodate annot        
        for j=1:nTrials
            disp(strcat(SUBJECT,' . ',chan_selected{e},' : bycycle + updating annoot tables for trial n:',string(j),'/',string(nTrials)))
            if isempty(pow_beta{j})
                continue
            end
            % 1) get bursts by freq domain general 
            cfg=[];
            cfg.threshold=percFooof;
            cfg.fs=TF_RATE;
            cfg.domain='frequency';
            maskTh=bursts_consistence(cfg,pow_beta{j});
            mask_beta{j}=maskTh;
            pieces=bwconncomp(maskTh);
            for ii=1:pieces.NumObjects
                idx=pieces.PixelIdxList{ii};
                
                % annot
                warning off
                annot_beta.id(bb)=ii;
                annot_beta.starts(bb)=Esingle_norm.time(idx(1));
                annot_beta.ends(bb)=Esingle_norm.time(idx(end));
                annot_beta.duration(bb)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1));
                annot_beta.trial_id(bb)=j;
                annot_beta.chan_id(bb)=chan_selected(e);
                annot_beta.power(bb)=1/TF_RATE*trapz(pow_beta{j}(idx)); 
                annot_beta.pow_max(bb)=max(pow_beta{j}(idx));
                annot_beta.pow_mean(bb)=mean(pow_beta{j}(idx),'omitnan');
                bb=bb+1;
                warning on
            end
            
            % 2) get bursts by freq domain iFreq from fooof
            cfg=[];
            cfg.threshold=percFooof;
            cfg.fs=TF_RATE;
            cfg.domain='frequency';
            maskTh=bursts_consistence(cfg,pow_fooof{j});
            mask_fooof{j}=maskTh;
            pieces=bwconncomp(maskTh);
            for ii=1:pieces.NumObjects
                idx=pieces.PixelIdxList{ii};
                % annot
                warning off
                annot_fooof.id(ff)=ii;
                annot_fooof.starts(ff)=Esingle_norm.time(idx(1));
                annot_fooof.ends(ff)=Esingle_norm.time(idx(end));
                annot_fooof.duration(ff)=Esingle_norm.time(idx(end))-Esingle_norm.time(idx(1)); 
                annot_fooof.trial_id(ff)=j;
                annot_fooof.chan_id(ff)=chan_selected(e);
                annot_fooof.power(ff)=1/TF_RATE*trapz(pow_fooof{j}(idx)); 
                annot_fooof.pow_max(ff)=max(pow_fooof{j}(idx));
                annot_fooof.pow_mean(ff)=mean(pow_fooof{j}(idx),'omitnan');
                annot_fooof.ifreq(ff)=ifreqs(e);
                ff=ff+1;
                warning on
            end

            % 3) get bursts by time domain by bycycle
            sig=data0(j,:);
                
            bycycle_results=bycycle(sig,fsample,f_range_bycycle,settings_bycycle,center_extrema,burst_method,return_model);
            bycycle_results.id=(1:height(bycycle_results))';
            bycycle_results.starts=Etime(bycycle_results.(first_sample) + 1)'; 
            bycycle_results.ends=Etime(bycycle_results.(last_sample) + 1)';            
            
            [labeledX, numRegions] = bwlabel(bycycle_results.is_burst);
            if numRegions>0 % is there burst?
                warning off
                for ii=1:numRegions
                    annot_cycle.id(tt)=tt;
                    annot_cycle.starts(tt)=min(bycycle_results.starts(labeledX == ii));
                    annot_cycle.ends(tt)=max(bycycle_results.ends(labeledX == ii));
                    annot_cycle.duration(tt)=max(bycycle_results.ends(labeledX == ii)) - min(bycycle_results.starts(labeledX == ii));
                    annot_cycle.trial_id(tt)=j;
                    annot_cycle.chan_id(tt)=chan_selected(e);
                    annot_cycle.volt_amp(tt)=mean(bycycle_results.volt_amp(labeledX == ii)); 
                    annot_cycle.volt_peak(tt)=mean(bycycle_results.volt_peak(labeledX == ii));
                    annot_cycle.band_amp(tt)=mean(bycycle_results.band_amp(labeledX == ii));
                    annot_cycle.time_rdsym(tt)=mean(bycycle_results.time_rdsym(labeledX == ii));
                    annot_cycle.time_ptsym(tt)=mean(bycycle_results.time_ptsym(labeledX == ii));
                    annot_cycle.cycles(tt)=sum(labeledX == ii);
                    annot_cycle.frequency(tt)=mean(bycycle_results.frequency(labeledX == ii));
                    tt=tt+1;  
                end 
                warning on
           end
        end
     
    end

%     %% add the window tag
%     disp('--adding tag window')
%     for ii=1:height(annot_beta)
%         burst_start=annot_beta.starts(ii);
%         burst_end=annot_beta.ends(ii);
%         burst_duration=burst_end-burst_start;
%         if (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end<=mean(tf_baseline_t0(:,2),'omitnan'))... % 1) area inside the window
%                 || (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end-burst_duration/2<=mean(tf_baseline_t0(:,2),'omitnan')) ... % 2) start inside and end outside only for half
%                 || (burst_start+burst_duration/2>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end<=mean(tf_baseline_t0(:,2),'omitnan')) ... % 3) start outside only for half and end inside
%                 || (burst_start<=mean(tf_baseline_t0(:,1),'omitnan') && burst_end>=mean(tf_baseline_t0(:,2),'omitnan')) % 4) window inside area
%             annot_beta.window(ii)={'baseline'}; 
%         elseif (burst_start>=x_stim1 && burst_end<=x_stim3_end) || (burst_start>=x_stim1 && burst_end-burst_duration/2<=x_stim3_end) || (burst_start+burst_duration/2>=x_stim1 && burst_end<=x_stim3_end) || (burst_start<=x_stim1 && burst_end>=x_stim3_end)
%             annot_beta.window(ii)={'stimulus'};
%         elseif (burst_start>=x_stim3_end && burst_end<=x_syl1) || (burst_start>=x_stim3_end && burst_end-burst_duration/2<=x_syl1) || (burst_start+burst_duration/2>=x_stim3_end && burst_end<=x_syl1) || (burst_start<=x_stim3_end && burst_end>=x_syl1)
%             annot_beta.window(ii)={'prespeech'};
%         elseif (burst_start>=x_syl1 && burst_end<=x_syl3_end) || (burst_start>=x_syl1 && burst_end-burst_duration/2<=x_syl3_end) || (burst_start+burst_duration/2>=x_syl1 && burst_end<=x_syl3_end) || (burst_start<=x_syl1 && burst_end>=x_syl3_end)
%             annot_beta.window(ii)={'speech'};
%         elseif (burst_start>=tf_rebound_t0(1) && burst_end<=tf_rebound_t0(2)) || (burst_start>=tf_rebound_t0(1) && burst_end-burst_duration/2<=tf_rebound_t0(2)) || (burst_start+burst_duration/2>=tf_rebound_t0(1) && burst_end<=tf_rebound_t0(2)) || (burst_start<=tf_rebound_t0(1) && burst_end>=tf_rebound_t0(2))
%             annot_beta.window(ii)={'rebound'};
%         else;annot_beta.window(ii)={'outside'};
%         end
%     end
% 
%     for ii=1:height(annot_fooof)
%         burst_start=annot_fooof.starts(ii);
%         burst_end=annot_fooof.ends(ii);
%         burst_duration=burst_end-burst_start;
%         if (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end<=mean(tf_baseline_t0(:,2),'omitnan'))... % 1) area inside the window
%                 || (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end-burst_duration/2<=mean(tf_baseline_t0(:,2),'omitnan')) ... % 2) start inside and end outside only for half
%                 || (burst_start+burst_duration/2>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end<=mean(tf_baseline_t0(:,2),'omitnan')) ... % 3) start outside only for half and end inside
%                 || (burst_start<=mean(tf_baseline_t0(:,1),'omitnan') && burst_end>=mean(tf_baseline_t0(:,2),'omitnan')) % 4) window inside area
%             annot_fooof.window(ii)={'baseline'}; 
%         elseif (burst_start>=x_stim1 && burst_end<=x_stim3_end) || (burst_start>=x_stim1 && burst_end-burst_duration/2<=x_stim3_end) || (burst_start+burst_duration/2>=x_stim1 && burst_end<=x_stim3_end) || (burst_start<=x_stim1 && burst_end>=x_stim3_end)
%             annot_fooof.window(ii)={'stimulus'};
%         elseif (burst_start>=x_stim3_end && burst_end<=x_syl1) || (burst_start>=x_stim3_end && burst_end-burst_duration/2<=x_syl1) || (burst_start+burst_duration/2>=x_stim3_end && burst_end<=x_syl1) || (burst_start<=x_stim3_end && burst_end>=x_syl1)
%             annot_fooof.window(ii)={'prespeech'};
%         elseif (burst_start>=x_syl1 && burst_end<=x_syl3_end) || (burst_start>=x_syl1 && burst_end-burst_duration/2<=x_syl3_end) || (burst_start+burst_duration/2>=x_syl1 && burst_end<=x_syl3_end) || (burst_start<=x_syl1 && burst_end>=x_syl3_end)
%             annot_fooof.window(ii)={'speech'};
%         elseif (burst_start>=tf_rebound_t0(1) && burst_end<=tf_rebound_t0(2)) || (burst_start>=tf_rebound_t0(1) && burst_end-burst_duration/2<=tf_rebound_t0(2)) || (burst_start+burst_duration/2>=tf_rebound_t0(1) && burst_end<=tf_rebound_t0(2)) || (burst_start<=tf_rebound_t0(1) && burst_end>=tf_rebound_t0(2))
%             annot_fooof.window(ii)={'rebound'};
%         else;annot_fooof.window(ii)={'outside'};
%         end
%     end
% 
%     for ii=1:height(annot_cycle)
%         burst_start=annot_cycle.starts(ii);
%         burst_end=annot_cycle.ends(ii);
%         burst_duration=burst_end-burst_start;
%         if (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end<=mean(tf_baseline_t0(:,2),'omitnan'))... % 1) area inside the window
%                 || (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end-burst_duration/2<=mean(tf_baseline_t0(:,2),'omitnan')) ... % 2) start inside and end outside only for half
%                 || (burst_start+burst_duration/2>=mean(tf_baseline_t0(:,1),'omitnan') && burst_end<=mean(tf_baseline_t0(:,2),'omitnan')) ... % 3) start outside only for half and end inside
%                 || (burst_start<=mean(tf_baseline_t0(:,1),'omitnan') && burst_end>=mean(tf_baseline_t0(:,2),'omitnan')) % 4) window inside area
%             annot_cycle.window(ii)={'baseline'}; 
%         elseif (burst_start>=x_stim1 && burst_end<=x_stim3_end) || (burst_start>=x_stim1 && burst_end-burst_duration/2<=x_stim3_end) || (burst_start+burst_duration/2>=x_stim1 && burst_end<=x_stim3_end) || (burst_start<=x_stim1 && burst_end>=x_stim3_end)
%             annot_cycle.window(ii)={'stimulus'};
%         elseif (burst_start>=x_syl1 && burst_end<=x_syl3_end) || (burst_start>=x_syl1 && burst_end-burst_duration/2<=x_syl3_end) || (burst_start+burst_duration/2>=x_syl1 && burst_end<=x_syl3_end) || (burst_start<=x_syl1 && burst_end>=x_syl3_end)
%             annot_cycle.window(ii)={'prespeech'};
%         elseif (burst_start>=x_stim3_end && burst_end<=x_syl1) || (burst_start>=x_stim3_end && burst_end-burst_duration/2<=x_syl1) || (burst_start+burst_duration/2>=x_stim3_end && burst_end<=x_syl1) || (burst_start<=x_stim3_end && burst_end>=x_syl1)
%             annot_cycle.window(ii)={'speech'};
%         elseif (burst_start>=tf_rebound_t0(1) && burst_end<=tf_rebound_t0(2)) || (burst_start>=tf_rebound_t0(1) && burst_end-burst_duration/2<=tf_rebound_t0(2)) || (burst_start+burst_duration/2>=tf_rebound_t0(1) && burst_end<=tf_rebound_t0(2)) || (burst_start<=tf_rebound_t0(1) && burst_end>=tf_rebound_t0(2))
%             annot_cycle.window(ii)={'rebound'};
%         else;annot_cycle.window(ii)={'outside'};
%         end
%     end

    %% save
    % correct id and save
    disp('-- saving')
    annot_beta.id(:)=1:height(annot_beta);filename=strcat(SUBJECT,'_beta_bursts.txt');
    writetable(annot_beta, strcat('annot/',filename), 'Delimiter', '\t');

    annot_fooof.id(:)=1:height(annot_fooof);filename=strcat(SUBJECT,'_fooof_bursts.txt');
    writetable(annot_fooof, strcat('annot/',filename), 'Delimiter', '\t');

    annot_cycle.id(:)=1:height(annot_cycle);filename=strcat(SUBJECT,'_cycle_bursts.txt');
    writetable(annot_cycle, strcat('annot/',filename), 'Delimiter', '\t');
end
end