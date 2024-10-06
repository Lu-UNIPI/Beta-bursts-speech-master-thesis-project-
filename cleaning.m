close all;clear all;clc;
%% 

PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long
% from previous analysis:
A=2; % median distance stim1 onset - syl 1 onset
B=1.26; % median distance syl1 onset - syl3 offset
IT=1.93; % median inter trial distance (syl3 offset - stim 1 onset)
LEFT=1.2;
RIGHT=1.5; 
visualization=0;

tStart=tic;
% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for i=3:3%numel(SUBJECTS) 
    % open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    % paths and load
    NAME_RAW_SIGNAL=strcat(SUBJECT,'_ft_raw_session');
    PATH_SIGNAL=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\FieldTrip\');
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');

    load(strcat(PATH_SIGNAL,filesep,NAME_RAW_SIGNAL));

    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    el_ecog=electrode(electrode.type=="ecog", {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1'});

    
    % take sessions with DBS LFP data
    id_session=loaded_epoch.session_id(strcmpi(loaded_epoch.type, 'LEAD'));

    %% data preprocessing
    % select channel and sessions - every session lasts around 5 min
    cfg=[];
    cfg.channel     ={'ecog*','dbs*'};
    cfg.trials      =id_session;
    D=ft_selectdata(cfg,D);
    D_preproc=preprocessing(D);

    %% choose channels to remove
    % bml_databrowser2
    cfg = [];
    cfg.channel = {'all','-ecog_115','-ecog_116','-ecog_124','-ecog_132','-ecog_152','-ecog_159','-ecog_163'}; 
    D_chan = ft_selectdata(cfg, D_preproc);

    %% referencing
    cfg=[];
    cfg.label   = el_ecog.electrode;
    cfg.group   = el_ecog.connector;
    cfg.method  = 'CTAR'; % using trimmed average referencing
    cfg.percent = 50; % percentage of 'extreme' channels in group to trim
    D_ref=bml_rereference(cfg,D_preproc);
    
    %% epoch
    % get timing info
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding'));
    coding=coding(coding.session_id==id_session, {'id', 'starts','ends','duration','session_id','trial_id','syl1_onset','syl3_offset'});
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise'));
    cue=cue(cue.session_id==id_session, {'id', 'starts','ends','duration','session_id','trial_id','stim1_starts'});

    % define trials
    % lock stimulus around syl1 onset lock=1, around stim1 onset lock=0
    epoch=table();
    epoch.id=(1:height(coding))';
    cfg=[];
    lock=true;
    if lock
        epoch.starts=coding.syl1_onset - (A+LEFT);
        epoch.ends=coding.syl1_onset + (B+RIGHT);
        epoch.syl1onset=coding.syl1_onset;
        % cfg.timelock='syl1onset';
    elseif ~lock
        epoch.starts=cue.stim1_starts - (IT+LEFT);
        epoch.ends=cue.stim1_starts + (A+RIGHT);
        epoch.stim1onset=cue.stim1_starts;
        % cfg.timelock='stim1onset';
    end
    cfg.epoch=epoch;
    [D_ecog_ep_art,info_epoch]=bml_redefinetrial(cfg,D_ref);

    [E_ref,info_epoch]=bml_redefinetrial(cfg,D_ref);
    save(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_ecog_trial_ref_globt.mat"), 'E_ref') 

    [E_preproc,info_epoch]=bml_redefinetrial(cfg,D_preproc);
    save(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_ecog_trial_globt.mat"), 'E_preproc') 

    %% choose trials to reject
    % bml_databrowser2

end