close all;clear all;clc;
%%
tab=readtable("HCPMMP1toAreas.txt");
area=struct();
for a=1:numel(unique(tab.area))
    label_area=unique(tab.area);
    label_area=label_area{a};
    area.(label_area)=0;
end
area.dbs=0;

%%
PATH_DATA='Z:\DBS';
% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

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
  
    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));

    % load data
    load(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_syl.mat"));
    %load(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_clean_stim.mat"));a=-(IT+LEFT);b=A+RIGHT;center=cue.stim1_starts(tokeepCue);
    

    %% select ecog and referencing
    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);
    
    % ecog
    choice='ecog';%'ecog','area'
    switch choice
        case 'area'
            channels=electrode.HCPMMP1_area;
            channels=channels(~(cellfun('isempty',channels)));
            idx=find(cellfun(@(x) strcmp(x,'AAC'),channels));
            % "PMC","POC","IFC","SMC","TPOJ","LTC","DLPFC","IFOC","AAC","IPC"
            if sum(idx)==0
                warning('No channels in this area');
            else
%                 chan_selected=E.label(idx);
                cfg=[];
                cfg.channel     =chan_selected;
                ECOG=ft_selectdata(cfg,E);
            end
        case 'ecog'
            cfg=[];
            cfg.channel     ={'ecog*'};
            ECOG=ft_selectdata(cfg,E);
    end

    cfg=[];
    cfg.label   = electrode.electrode;
    cfg.group   = electrode.connector;
    cfg.method  = 'CTAR';   % using trimmed average referencing
    cfg.percent = 50;       % percentage of 'extreme' channels in group to trim
    ECOG=bml_rereference(cfg,ECOG);
    % save(fullfile(PATH_SIGNAL,'LFP beta burst',SUBJECT+"_rebound_clean_ref_globt.mat"), 'E_ctar')
   

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
        if startsWith(chan_selected{e},'dbs');area.dbs=area.dbs+1;continue;end
        label_area=string(electrode.HCPMMP1_area(strcmp(electrode.electrode,chan_selected{e})));
        label_area=label_area(1);
        if strcmp(label_area,'out of atlas') || strcmp(label_area,'label not present');continue;end
        area.(label_area)=area.(label_area)+1;
    end
end

%%
for a=1:(numel(unique(tab.area))+1)   
    if a==1
        area_label='dbs';
        color=[1,0.95,1];
    else
        area_short=unique()
        area_label=tab.area_fullname(strcmp(tab.area,tab.area(a)));area_label=string(area_label(1));
        color=tab.color(strcmp(tab.area,tab.area(a)));color=string(color(1));
    end
    histogram()
end
histogram()