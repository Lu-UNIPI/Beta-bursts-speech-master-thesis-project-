clear all;close all;clc
%%
% This script finds median in timing (A and B) to define epochs
% Epoch lock in in speech onset will be inside A - B + borders
% Epoch lock in in stimulus presentation will be inside IT - A + borders

PATH_DATA='Z:\DBS';
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

disp(repelem('*',100))
diary log_epoch_timing
disp("A =: time window between stim1 onset and syl1 onset")
disp("B =: time window between syl1 onset and syl3 offset")

tStart=tic;

epoch=[];
epoch.A=cell(1,numel(SUBJECTS)); % time window stim1 - syl1
epoch.B=cell(1,numel(SUBJECTS)); % time window syl1 - syl3
for i=1:numel(SUBJECTS)
    % open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    sessions=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    
    % take sessions with DBS LFP data
    id_session=sessions.id(strcmpi(sessions.type, 'LEAD'));

    % take start and stopping trial info from annot tables
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding'));
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise'));
    timing_cue=cue(cue.session_id==id_session,{'trial_id','stim1_starts'});
    timing_coding=coding(coding.session_id==id_session,{'trial_id','syl1_onset','syl3_offset'});
    
    % join correctly tables
    timing_trials=join_timing(timing_cue,timing_coding);
    
    % find duration breaks bewtween trials
    timing_trials.A=timing_trials.syl1_onset - timing_trials.stim1_starts;
    timing_trials.B=timing_trials.syl3_offset - timing_trials.syl1_onset;

    disp(SUBJECT+ '  median A: '+string( median(timing_trials.A,'omitnan')) )
    disp(SUBJECT+ '  median B: '+string( median(timing_trials.B,'omitnan')) )
    disp(repelem('-',25))
    epoch.A{1,i}=timing_trials.A;
    epoch.B{1,i}=timing_trials.B;

end
tEnd=toc(tStart);

%% distribution 
toplot_A = cellfun(@(x) x(:)', epoch.A, 'UniformOutput', false);
toplot_A = [toplot_A{:}];
toplot_B = cellfun(@(x) x(:)', epoch.B, 'UniformOutput', false);
toplot_B = [toplot_B{:}];
histogram(toplot_A);hold on;histogram(toplot_B)
text(0.7*max(xlim),0.8*max(ylim),('Median A: '+string(median(toplot_A,'omitnan'))))
text(0.7*max(xlim),0.6*max(ylim),('Median B: '+string(median(toplot_B,'omitnan'))))
title('Distribution of timing before and after syl1 onset for PD with DBS data')
saveas(gcf,'Timing around syl1 onset.png')
disp('Median tot A:   '+string(median(toplot_A,'omitnan')))
disp('Median tot B:   '+string(median(toplot_B,'omitnan')))

diary off
%% definitions
% BASELINE -1s stim1_start
% REBOUND -300ms offset_syl3 +700ms


