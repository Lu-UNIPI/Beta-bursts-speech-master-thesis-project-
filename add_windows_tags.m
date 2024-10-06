close all; clear all; clc

%% add the window tags
PATH_DATA='Z:\DBS';
PATH_ANNOT_internal = 'annot/general CTAR/';
files = dir([pwd,'/annot/general CTAR/*.txt']);
for file = files'
    tab = readtable([PATH_ANNOT_internal,file.name]);

    if ismember('window50', tab.Properties.VariableNames);continue;end

    SUBJECT=file.name(1:7);
    disp(['Now running ',file.name])
    
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

    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);

    tokeepCod=(~isnan(coding.syl1_onset));
    tokeepCue=(~isnan(cue.stim1_starts));
    switch SUBJECT(4:end)
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

    windows={'baseline','stimulus','prespeech','speech','rebound'};
    timing_windows=struct();
    for w=1:numel(windows)
        window=windows{w};
        switch window
            case 'baseline'
                timing_windows.(window).start= mean(tf_baseline_t0(:,1),'omitnan');
                timing_windows.(window).end= mean(tf_baseline_t0(:,2),'omitnan');
            case 'stimulus'
                timing_windows.(window).start= x_stim1;
                timing_windows.(window).end= x_stim3_end;
            case 'prespeech'
                timing_windows.(window).start= x_stim3_end;
                timing_windows.(window).end= x_syl1;
            case 'speech'
                timing_windows.(window).start= x_syl1;
                timing_windows.(window).end= x_syl3_end;
            case 'rebound'
                timing_windows.(window).start= tf_rebound_t0(1);
                timing_windows.(window).end= tf_rebound_t0(2);
        end
    end

    %% change window tag already saved
    disp("--removing existing window's tag ")
    if ismember('window', tab.Properties.VariableNames)
        tab=removevars(tab, 'window');
    end

    % add other possible window tags
    %% 60,50,40,30,20
    window_tags = {'window60','window50','window40','window30','window20'};
    
    for t=1:numel(window_tags)
        tag=window_tags{t};
        disp(['--adding tag ',tag])
        percentage= str2double(tag(end-1:end));
        for ii=1:height(tab)
            burst_start=tab.starts(ii);
            burst_end=tab.ends(ii);
            burst_duration=burst_end-burst_start;

            for w = 1:numel(windows)
                window=windows{w};
                A=timing_windows.(window).start;
                B=timing_windows.(window).end;
                if ( burst_start >= A && burst_end <= B)... % 1) area inside the window
                        || (burst_start >= A && burst_end - burst_duration*(100-percentage)/100 <= B) ... % 2) starts inside and only percentage % is inside, and ends outside
                        || (burst_start + burst_duration*(100-percentage)/100 >= A && burst_end <= B) ... % 3) starts outside but percentage % is inside, and ends inside
                        || (burst_start <= A && burst_end >= B) % 4) window inside area
                    tab.(tag)(ii)={window};
                    
                    % check if the burst is present also in the next window. If yes, decide on to which belong the majority of the burst
                    if ~(w==numel(windows))
                        % percentage of burst inside the window
                        if burst_start < A; perc_start=A; else; perc_start=burst_start;end
                        if burst_end > B; perc_end=B; else; perc_end=burst_end;end
                        perc_in = (perc_end - perc_start) / burst_duration;
                        
                        C=timing_windows.(windows{w+1}).start;
                        D=timing_windows.(windows{w+1}).end;
                        if ( burst_start >= C && burst_end <= D) || (burst_start >= C && burst_end - burst_duration*(100-percentage)/100 <= D) || (burst_start + burst_duration*(100-percentage)/100 >= C && burst_end <= D) || (burst_start <= C && burst_end >= D) 
                            % percentage of burst inside the window
                            if burst_start < C; perc_start=C; else; perc_start=burst_start;end
                            if burst_end > D; perc_end=D; else; perc_end=burst_end;end
                            perc_in_next = (perc_end - perc_start) / burst_duration;
                            if perc_in_next > perc_in; tab.(tag)(ii)={windows{w+1}};end
                        end
                    end
                end
            end  
        end
        tab.(tag)(cellfun(@isempty, tab.(tag))) = {'outside'};
    end


    %% onset
    disp('--adding tag onset')
    tag='window_onset';
    for ii=1:height(tab)
        burst_start=tab.starts(ii);
        burst_end=tab.ends(ii);
        burst_duration=burst_end-burst_start;
        if (burst_start>=mean(tf_baseline_t0(:,1),'omitnan') && burst_start<=mean(tf_baseline_t0(:,2),'omitnan')) % start is inside the window
            tab.(tag)(ii)={'baseline'}; 
        elseif (burst_start>=x_stim1 && burst_start<=x_stim3_end) 
            tab.(tag)(ii)={'stimulus'};
        elseif (burst_start>=x_stim3_end && burst_start<=x_syl1) 
            tab.(tag)(ii)={'prespeech'};
        elseif (burst_start>=x_syl1 && burst_start<=x_syl3_end) 
            tab.(tag)(ii)={'speech'};
        elseif (burst_start>=tf_rebound_t0(1) && burst_start<=tf_rebound_t0(2)) 
            tab.(tag)(ii)={'rebound'};
        else;tab.(tag)(ii)={'outside'};
        end
    end

    writetable(tab,[PATH_ANNOT_internal,file.name], 'Delimiter', '\t');
end
disp('Finish :)')