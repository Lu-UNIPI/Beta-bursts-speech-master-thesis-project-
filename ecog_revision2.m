close all;clear all;clc;
%% call struct and open
% here you can have the original data, the filtered ones, and the rejected ones
% Your consideration can be saved in a diary log file

PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

for i=2:2%numel(SUBJECTS) 
    % paths and load
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    NAME_DATA_FILT=strcat(SUBJECT,'_ecog_trial_critA_ref_globt.mat'); 
    NAME_DATA_ORIGINAL=strcat(SUBJECT,'_ecog_trial_globt.mat');
    PATH_SIGNAL=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\FieldTrip\LFP beta burst\');
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    load(fullfile(PATH_SIGNAL,NAME_DATA_FILT));
    load(fullfile(PATH_SIGNAL,NAME_DATA_ORIGINAL));


    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    
    %electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    %el_ecog=electrode(electrode.type=="ecog", {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1'});
    
    artifacts_critA=bml_annot_read(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_ecog_trial_critA_ref_globt.txt')));
    artifacts_critA=artifacts_critA(:,1:end-1);
    info_epoch=bml_annot_read(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_info_epoch_syl1.txt')));
    continue
    %% Manual artifact rejection
    % bml_databrowser2

    % remove single channels single trials
    try
        tab_torestore=readtable(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_torestore.txt')));
        artifacts_torestore=tab_to_annot(tab_torestore,info_epoch);
        bml_annot_write(artifacts_torestore,fullfile(PATH_ANNOT,'LFP beta burst',SUBJECT+"_ecog_trial_critA_ref_torestore_globt"))

        tab_toreject=readtable(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_toreject.txt')));
        artifacts_toreject=tab_to_annot(tab_toreject,info_epoch);
        bml_annot_write(artifacts_toreject,fullfile(PATH_ANNOT,'LFP beta burst',SUBJECT+"_ecog_trial_critA_ref_toreject_globt"))
        
        cfg=[];
        cfg.groupby='label';
        a=bml_annot_difference(cfg,artifacts_critA,artifacts_torestore);
        artifacts_final=sortrows(bml_annot_rowbind(a,artifacts_toreject),'starts');
        artifacts_final.id=(1:height(artifacts_final))';
        bml_annot_write(artifacts_final,fullfile(PATH_ANNOT,'LFP beta burst',SUBJECT+"_ecog_trial_critA_ref_manual_globt"))
        
        
    catch
        tab_manual=readtable(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_manual.txt')));
        artifacts_final=tab_to_annot(tab_manual,info_epoch);
        bml_annot_write(artifacts_final,fullfile(PATH_ANNOT,'LFP beta burst',SUBJECT+"_ecog_trial_critA_ref_manual_globt"))

    end

    cfg=[];
    cfg.label_colname   = 'label';	        % column of annot with label specification
    cfg.annot           = artifacts_final;  % annotation table to use for mask
    cfg.value           = NaN;              % value to use for mask
    E=bml_mask(cfg, Ebase);
    E=consolidate(E,Ebase);

    % remove full channels
    try 
        chans=readmatrix(fullfile(PATH_ANNOT,'LFP beta burst',strcat(SUBJECT,'_channelstoremove.txt')));
        chans=arrayfun(@(x) ['ecog_',num2str(x)],chans,'UniformOutput', false);
        labels=setdiff(unique(E.label),chans);
        cfg=[];
        cfg.channel     =labels;
        E=ft_selectdata(cfg,E);
    catch
        disp('No chans to remove')
    end

    save(fullfile(PATH_SIGNAL,SUBJECT+"_ecog_trial_critA_ref_manual_globt.mat"), 'E')
    
end
