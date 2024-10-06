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
delta=1.1;

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
settings_bycycle = struct();
settings_bycycle.amp_fraction=0.4;
settings_bycycle.amp_consistency=0.35;
settings_bycycle.period_consistency=0.55;
settings_bycycle.monotonicity=0.6;
f_range_bycycle = [10,40];
return_model = false;

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for cc=1
for i=5%:numel(SUBJECTS) 
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
    prod_triplet=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_produced_triplet'));prod_triplet=prod_triplet(prod_triplet.session_id==id_session,:);
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
    % annot_beta=table();
    
    % remasking NaNs with zeros
    cfg=[];
    cfg.value           =0;
    cfg.remask_nan      =true;
    cfg.complete_trial  =true;
    E=bml_mask(cfg,E);

    nElec=numel(chan_selected); 
    ifreqs=nan(1,nElec);
    for e=31:nElec
        disp('---------------------------------------------------------')
        fprintf('%s %s',SUBJECT,chan_selected{e})
        pow_beta={};mask_beta={};
        pow_fooof={};mask_fooof={};
        pow_time={};mask_time={};mask_time_toplot={};
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
            if sum(data(:)==0)/numel(data)>=0.7
                pow_beta{j}=[];pow_fooof{j}=[];
                %if j==nTrials;pow_beta{j}=[];pow_fooof{j}=[];ifreqs(j)=NaN;end
                disp('Trial rejected, going to the next');
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

            % 1) get poi in domain general
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

        %% FIGURE: plot signals 
        fig1=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
        subplot(3,3,[1,4])
        xlabel('Time [s]')
        ylabel('Power [V^2]')
        title('Frequency based - general beta')
        hold on
        subplot(3,3,[2,5])
        xlabel('Time [s]')
        ylabel('Power [V^2]')
        title('Frequency based - fooof iFreq')
        hold on
        subplot(3,3,[3,6])
        xlabel('Time [s]')
        ylabel('Power [V^2]')
        title('Time based - byclycle')
        hold on
        sgtitle({[SUBJECT+": "+strrep(chan_selected{e},'_',' ')+' '+'bursts identification methods comparison'],''})
        for j=1:nTrials
            if isempty(pow_beta{j})
                mask_beta{j}=logical(zeros(size(Esingle_norm.time)));
                mask_fooof{j}=logical(zeros(size(Esingle_norm.time)));
                mask_time{j}=logical(zeros(size(Etime)));
                mask_time_toplot{j}=logical(zeros(size(Esingle_norm.time)));
                continue
            end
            disp(strcat(SUBJECT,' . ',chan_selected{e},' : bycycle + preparing plot for trial n:',string(j),'/',string(nTrials)))

            % 1) get bursts by freq domain general
            subplot(3,3,[1,4])
            ypos=delta*(j-1); 
            plot(Esingle_norm.time,pow_beta{j}+ypos,'b');hold on
            %text(Esingle_norm.time(1)-0.3, ypos , num2str(j)); 
            cfg=[];
            cfg.threshold=percBeta;
            cfg.fs=TF_RATE;
            cfg.domain='frequency';
            maskTh=bursts_consistence(cfg,pow_beta{j});
            mask_beta{j}=maskTh;
            pieces=bwconncomp(maskTh);
            for ii=1:pieces.NumObjects
                idx=pieces.PixelIdxList{ii};
                plot(Esingle_norm.time(idx),pow_beta{j}(idx)+ypos,'r','Linewidth',1.5);hold on
            end
            xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1.5)
            xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1.5)
            xline([x_stim1 x_stim1],'Color','g','LineWidth',1.5)
            %xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
            %xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
            xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1.5)
            xline([x_syl1 x_syl1],'Color','r','LineWidth',1.5)
            %xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
            %xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
            xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1.5)
            ymax=max(pow_beta{j}+ypos);
            ylim([-1 ymax]);
            % xlim([Esingle_norm.time(1)-0.5,Esingle_norm.time(end)])
            xlim([Esingle_norm.time(1),Esingle_norm.time(end)])

            % 2) get bursts by freq domain iFreq from fooof
            subplot(3,3,[2,5]) 
            ypos=delta*(j-1);
            plot(Esingle_norm.time,pow_fooof{j}+ypos,'b');hold on
            %text(Esingle_norm.time(1)-0.3, ypos , num2str(j)); 
            cfg=[];
            cfg.threshold=percFooof;
            cfg.fs=TF_RATE;
            cfg.domain='frequency';
            maskTh=bursts_consistence(cfg,pow_fooof{j});
            mask_fooof{j}=maskTh;
            pieces=bwconncomp(maskTh);
            for ii=1:pieces.NumObjects
                idx=pieces.PixelIdxList{ii};
                plot(Esingle_norm.time(idx),pow_fooof{j}(idx)+ypos,'r','Linewidth',1.5);hold on
            end
            xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1.5)
            xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1.5)
            xline([x_stim1 x_stim1],'Color','g','LineWidth',1.5)
            %xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
            %xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
            xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1.5)
            xline([x_syl1 x_syl1],'Color','r','LineWidth',1.5)
            %xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
            %xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
            xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1.5)
            ymax=max(pow_fooof{j}+ypos);
            ylim([-1 ymax]);
            xlim([Esingle_norm.time(1),Esingle_norm.time(end)])

            % 3) get bursts by time domain by bycycle    
            sig=data0(j,:);
            bycycle_results=bycycle(sig,fsample,f_range_bycycle,settings_bycycle,center_extrema,burst_method,return_model);
            bycycle_results.id=(1:height(bycycle_results))';
            bycycle_results.starts=Etime(bycycle_results.(first_sample) + 1)'; 
            bycycle_results.ends=Etime(bycycle_results.(last_sample) + 1)'; 
            [labeledX, numRegions] = bwlabel(bycycle_results.is_burst);
            
            subplot(3,3,[3,6])
            
%             figure;
%             plot(Etime,sig+ypos,'b','LineWidth',0.1);hold on
%             maskBy=zeros(size(Etime));
%             for ii=1:numRegions
%                 time_start=min(bycycle_results.starts(labeledX == ii));
%                 [~,idx_start]=min(abs(Etime-time_start));
%                 time_end=max(bycycle_results.ends(labeledX == ii));
%                 [~,idx_end]=min(abs(Etime-time_end));
%                 plot(Etime(idx_start:idx_end),sig(idx_start:idx_end)+ypos,'r','Linewidth',1.5);hold on
% 
%                 maskBy(idx_start:idx_end)=1;
%             end



            ypos=delta*(j-1);
            plot(Esingle_norm.time,pow_beta{j}+ypos,'b','LineWidth',0.1);hold on
            %text(Esingle_norm.time(1)-0.3,ypos,num2str(j));
            maskBy=zeros(size(Esingle_norm.time));
            for ii=1:numRegions
                time_start=min(bycycle_results.starts(labeledX == ii));
                [~,idx_start]=min(abs(Esingle_norm.time-time_start));
                time_end=max(bycycle_results.ends(labeledX == ii));
                [~,idx_end]=min(abs(Esingle_norm.time-time_end));
                plot(Esingle_norm.time(idx_start:idx_end),pow_beta{j}(idx_start:idx_end)+ypos,'r','Linewidth',1.5);hold on

                maskBy(idx_start:idx_end)=1;
            end 
            xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1.5)
            xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1.5)
            xline([x_stim1 x_stim1],'Color','g','LineWidth',1.5)
            %xline([x_stim2 x_stim2],'Color','g','LineWidth',1,'LineStyle',':')
            %xline([x_stim3 x_stim3],'Color','g','LineWidth',1,'LineStyle',':')
            xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1.5)
            xline([x_syl1 x_syl1],'Color','r','LineWidth',1.5)
            %xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
            %xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
            xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1.5)
            ymax=max(pow_beta{j}+ypos);
            ylim([-1,ymax]) 
            xlim([Esingle_norm.time(1),Esingle_norm.time(end)])

            mask_time_toplot{j}=logical(maskBy);
       
        end
        
        % plot masks
        maskBeta=cell2mat(mask_beta');
        maskFooof=cell2mat(mask_fooof');
        maskTime=cell2mat(mask_time_toplot');
        for z=1:3
            subplot(3,3,6+z)
            switch z
                case 1;matrix=maskBeta;
                case 2;matrix=maskFooof;
                case 3;matrix=maskTime;
            end
            imagesc(Esingle_norm.time,1:nTrials,matrix)
            xline([mean(tf_baseline_t0(:,1),'omitnan') mean(tf_baseline_t0(:,1),'omitnan')],'Color','y','LineWidth',1.5)
            xline([mean(tf_baseline_t0(:,2),'omitnan') mean(tf_baseline_t0(:,2),'omitnan')],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(1) tf_rebound_t0(1)],'Color','y','LineWidth',1.5)
            xline([tf_rebound_t0(2) tf_rebound_t0(2)],'Color','y','LineWidth',1.5)
            xline([x_stim1 x_stim1],'Color','g','LineWidth',1.5)
            %xline([x_stim2 x_stim2],'Color','g','LineWidth',2,'LineStyle',':')
            %xline([x_stim3 x_stim3],'Color','g','LineWidth',2,'LineStyle',':')
            xline([x_stim3_end x_stim3_end],'Color','g','LineWidth',1.5)
            xline([x_syl1 x_syl1],'Color','r','LineWidth',1.5)
            %xline([x_syl2 x_syl2],'Color','r','LineWidth',1,'LineStyle',':')
            %xline([x_syl3 x_syl3],'Color','r','LineWidth',1,'LineStyle',':')
            xline([x_syl3_end x_syl3_end],'Color','r','LineWidth',1.5)
            axis square;
            xlabel('Time [s]')
            ylabel('Power [V^2]')
        end
        
        
        %% FIGURE: burst overlap
        maskBeta_line=reshape(maskBeta', [], 1)';
        maskFooof_line=reshape(maskFooof', [], 1)';
        maskTime_line=reshape(maskTime', [], 1)';
        len_mask=length(maskBeta_line);

        n_perm=10^3;
        p_values=nan(1,9);
        to_test=nan(1,9);
        fig2=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
        disp('-----------------------------------------')
        for z=1:3
            switch z
                %case 1;first=maskBeta_line;second=maskBeta_line;tit='Original beta vs original beta';
                %case 2;first=maskFooof_line;second=maskFooof_line;tit='Fooof vs fooof';
                %case 3;first=maskTime_line;second=maskTime_line;tit='Cycle vs cycle';
                case 1;first=maskBeta_line;second=maskFooof_line;tit='Original beta vs fooof';
                %case 5;first=maskFooof_line;second=maskBeta_line;tit='Fooof vs original beta';
                %case 6;first=maskTime_line;second=maskBeta_line;tit='Cycle vs original beta';
                case 2;first=maskBeta_line;second=maskTime_line;tit='Original beta vs cycle';
                case 3;first=maskFooof_line;second=maskTime_line;tit='Fooof vs cycle';
                %case 9;first=maskTime_line;second=maskFooof_line;tit='Cycle vs fooof';
            end
            disp(strcat(string(z),' : ',tit))
            cfg=[];
            cfg.n_permutations=n_perm;
            newmasks=create_mask(cfg,first);
            disp('Masks created')

            null_distribution=nan(1,n_perm);
            for perm=1:n_perm;null_distribution(perm)=length(intersect(find(first),find(newmasks(perm,:))))/len_mask;end
            disp('Null distribution created')

            subplot(1,3,z)
            histogram(null_distribution,7);hold on
            to_test(z)=length(intersect(find(first),find(second)))/len_mask;
            xline(to_test(z))
            title(tit)
            p_values(z)=sum(null_distribution>to_test(z))/n_perm;
            xlim([0, to_test(z)*1.1])
            xlabel('Overlap fraction')
            ylabel('Probability density')
        end
        sgtitle(strcat('Burst overlap comparison with ',string(n_perm),'permutations'))
        for z=1:9
            subplot(1,3,z);
            xlim([0, max(to_test)*1.1]); 
        end
        
        

        %% FIGURE: mask per condition
        % plot masks
        [~,idx_vol]=sort(prod_triplet.rms_audio_p(tokeepCod),'descend');
        [~,idx_rt]=sort(coding.syl3_offset(tokeepCod)-coding.syl1_onset(tokeepCod),'descend');
        try [~,idx_acc]=sort(prod_triplet.phonetic_ontarget_perc(tokeepCod),'descend');catch;end

        fig3=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
        for z=1:12
            subplot(4,3,z)
            switch z
                case 1;imagesc(Esingle_norm.time,1:nTrials,maskBeta);tit='Original beta';
                case 2;imagesc(Esingle_norm.time,1:nTrials,maskFooof);tit='Original fooof';
                case 3;imagesc(Esingle_norm.time,1:nTrials,maskTime);tit='Original cycle';
                case 4;imagesc(Esingle_norm.time,idx_vol,maskBeta(idx_vol,:));tit='Beta by volume';
                case 5;imagesc(Esingle_norm.time,idx_vol,maskFooof(idx_vol,:));tit='Fooof by volume';
                case 6;imagesc(Esingle_norm.time,idx_vol,maskTime(idx_vol,:));tit='Cycle by volume';
                case 7;imagesc(Esingle_norm.time,idx_rt,maskBeta(idx_rt,:));tit='Beta by RT';
                case 8;imagesc(Esingle_norm.time,idx_rt,maskFooof(idx_rt,:));tit='Fooof by RT';
                case 9;imagesc(Esingle_norm.time,idx_rt,maskTime(idx_rt,:));tit='Cycle by RT';
                case 10;try imagesc(Esingle_norm.time,idx_acc,maskBeta(idx_acc,:));tit='Beta by accuracy';catch;end
                case 11;try imagesc(Esingle_norm.time,idx_acc,maskFooof(idx_acc,:));tit='Fooof by accuracy';catch;end
                case 12;try imagesc(Esingle_norm.time,idx_acc,maskTime(idx_acc,:));tit='Cycle by accuracy';catch;end
            end
            
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
            axis square;
            title(tit)
            xlabel('Time [s]')
            ylabel('Power [V^2]')
        end
        sgtitle('Burst comparison per condition')

        %% 
        fig4=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
        plot(Esingle_norm.time,mean(maskBeta,1,'omitnan'));hold on
        plot(Esingle_norm.time,mean(maskFooof,1,'omitnan'))
        plot(Esingle_norm.time,mean(maskTime,1,'omitnan'))
        legend('General beta','Fooof iFreq','Bycycle')
        xlabel('Time [s]')
        ylabel('Burst probability')
        xlim([-3.4,2.7])
        title(strcat(strrep(string(Esingle_norm.label),'_','')," : ",'Burst probability over time'))

        %% saving
        name=[strcat('images/bursts CTAR/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_bursts_signals2')];
        saveas(fig1,strcat(name,'.png'))
        saveas(fig1,strcat(name,'.fig'))

        name=[strcat('images/bursts CTAR/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_bursts_overlap3')];
        saveas(fig2,strcat(name,'.png'))
        saveas(fig2,strcat(name,'.fig'))

        name=[strcat('images/bursts CTAR/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_bursts_condition')];
        saveas(fig3,strcat(name,'.png'))
        saveas(fig3,strcat(name,'.fig'))

        name=[strcat('images/bursts CTAR/',SUBJECT,'_speech_',strrep(string(Esingle_norm.label),'_',''),'_bursts_probability')];
        saveas(fig4,strcat(name,'.png'))
        saveas(fig4,strcat(name,'.fig'))

        close all

    end
end
end

%             bycycle_bursts=bycycle_results(bycycle_results.is_burst==1,:);
%             bycycle_bursts.id(:)=1:height(bycycle_bursts);
%             maskBy=false(1,length(Etime));
%             for z=1:height(bycycle_bursts);maskBy(bycycle_bursts.(first_sample)(z)+1: bycycle_bursts.(last_sample)(z)+1)=1;end
% 
%             cfg=[];
%             cfg.fs=fsample;
%             cfg.domain='time';
%             maskBy=bursts_consistence(cfg,maskBy);
%             mask_time{j}=maskBy;
%             maskBy_plot=false(1,length(pow_beta{1}));
%             for pp=1:length(maskBy_plot)
%                 [~,idx]=min(abs(Etime-Esingle_norm.time(pp)));
%                 maskBy_plot(pp)=maskBy(idx);
%             end
%             mask_time_toplot{j}=maskBy_plot;
%             [but_b,but_a]=butter(4, f_range/(fsample / 2), 'bandpass');
%             but_sig=filtfilt(but_b,but_a,sig);

%             sig=pow_beta{j};
%             maskBy=maskBy_plot;
%             deltat=delta;
%             Etime=Esingle_norm.time;


%             subplot(3,3,[3,6])
%             ypos=deltat*(j-1);
%             plot(Etime,sig+ypos,'b','LineWidth',0.1);hold on
%             text(Etime(1)-0.3,ypos,num2str(j));
%             pieces=bwconncomp(maskBy);
%             for ii=1:pieces.NumObjects
%                 idx=pieces.PixelIdxList{ii};
%                 plot(Etime(idx),sig(idx)+ypos,'r','Linewidth',1);hold on
%             end 
%             ymax=max(sig+ypos);
%             ylim([-100,ymax]) 
%             xlim([Etime(1)-0.5,Etime(end)])