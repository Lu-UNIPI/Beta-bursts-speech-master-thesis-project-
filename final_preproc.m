close all;clear all;clc;
%% 
%%% this script get raw data and 
% 1) select session id with DBS LFP data
% 2) preprocess 
% 3) epoch
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


% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for cc=1:2 % 1 for syl locked and 2 for stim locked
for i=14%1:numel(SUBJECTS) 
    %open a subject dir
    if i==6;continue;end
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    % paths and load
    NAME_RAW_SIGNAL=strcat(SUBJECT,'_ft_raw_session');
    PATH_SIGNAL=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\FieldTrip\');
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    disp('-- loading data')
    load(strcat(PATH_SIGNAL,filesep,NAME_RAW_SIGNAL));

    % take sessions with DBS LFP data
    id_session=loaded_epoch.session_id(strcmpi(loaded_epoch.type, 'LEAD'));

    %% preprocessing
    % select channel and sessions - every session lasts around 5 min
    disp('-- preprocessing')
    cfg=[];
    cfg.channel     ={'ecog*','dbs_L*'};
    cfg.trials      =id_session;
    D=ft_selectdata(cfg,D);
    D_preproc=preprocessing(D);

    %% correct coding and cue tables
    % get timing info
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding'));
    coding=coding(coding.session_id==id_session,:);
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise'));
    cue=cue(cue.session_id==id_session, :);

    tokeepCod=(~isnan(coding.syl1_onset));
    tokeepCue=(~isnan(cue.stim1_starts));
    switch string(SUBJECTS(i))
        case {'3003', '3010', '3015', '3020', '3025', '3027'};tokeepCue=tokeepCod;
        case {'3014','3012'};tokeepCod=tokeepCue;
        case {'3024'};tokeepCod=tokeepCod & tokeepCue;tokeepCue=tokeepCod;
        case {'3021'};idx=find(~tokeepCue)+1;idx=idx(1:2);tokeepCod(idx)=0;
    end

    %% epoching
    disp('-- epoching')
    epoch=table();
    epoch.id=(1:sum(tokeepCod))';
    cfg=[];
    switch cc
        case 1;center=coding.syl1_onset(tokeepCod);
        case 2;center=cue.stim1_starts(tokeepCue);
    end
    epoch.starts=center - (A+LEFT);
    epoch.ends=center + (B+RIGHT);
    epoch.syl1onset=center;
    epoch=bml_annot_extend(epoch,1);
    %cfg.timelock='syl1onset';
    cfg.epoch=epoch;
    [Ebase,info_epoch]=bml_redefinetrial(cfg,D_preproc);
    
    %% artifacts removal 
    disp('-- removing artifacts')
    if string(SUBJECTS(i))~="3012"
    if cc==1 
        tab_artifacts=readtable(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_manual_artifacts.txt')));
        artifacts=tab_to_annot(tab_artifacts,info_epoch);
        if exist(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_dbsmanual_artifacts.txt')),'file')
            tab_dbsartifacts=readtable(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_dbsmanual_artifacts.txt')));
            artifacts_dbs=tab_to_annot(tab_dbsartifacts,info_epoch);
            artifacts=[artifacts;artifacts_dbs];
        end
        artifacts=unique(artifacts(:,2:end),'rows');artifacts.id=(1:height(artifacts))';
        %artifacts=sortrows(artifacts,{'label','trial_id'});
        %artifacts.id=(1:height(artifacts))';
        bml_annot_write(artifacts,fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_artifacts.txt')));
    elseif cc==2
        artifacts=bml_annot_read(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_artifacts.txt')));
    end
    
    % remove artifacts
    cfg=[];
    cfg.label_colname   = 'label';	        % column of annot with label specification
    cfg.annot           = artifacts;        % annotation table to use for mask
    cfg.value           = NaN;              % value to use for mask
    E=bml_mask(cfg,Ebase);
    E=consolidate(E,Ebase);
    
     end
    
    %% SAVE 
    % convert in locked timing
    disp('-- saving')
    for t=1:length(center);E.time{t}=Ebase.time{t}-center(t);end
    if cc==1
        save(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_syl.mat"), 'E') 
    elseif cc==2
        save(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_stim.mat"), 'E') 
    end
    
end
end
disp('Finish! :)')

    
