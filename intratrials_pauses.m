clear all;close all;clc
%%
% This script define IT median time, useful to define epochs, baseline and
% rebound periods

PATH_DATA='Z:\DBS';
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

disp(repelem('*',100))
diary log_intratrial_pauses
disp("Displaying all medians: ")

tStart=tic;

durations=cell(1,numel(SUBJECTS));
n_trials=cell(1,numel(SUBJECTS));
for i=1:numel(SUBJECTS)
    % open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    sessions=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    
    % take sessions with DBS LFP data
    id_session=sessions.id(strcmpi(sessions.type, 'LEAD'));

    % take start and stopping trial info from annot tables
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding'));
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise'));
    timing_starts=cue(cue.session_id==id_session,{'trial_id','stim1_starts'});
    timing_ends=coding(coding.session_id==id_session,{'trial_id','syl3_offset'});
    
    % join correctly tables
    timing_trials=join_timing(timing_starts,timing_ends);
    
    % find duration breaks bewtween trials
    timing_trials.pause=[NaN; timing_trials.stim1_starts(2:end) - timing_trials.syl3_offset(1:end-1)];

    disp(SUBJECT+ ':  '+string( median(timing_trials.pause,'omitnan')) )
    durations{1,i}=timing_trials.pause;
    n_trials{1,i}=height(timing_trials);
end
tEnd=toc(tStart);

%% distribution 
toplot = cellfun(@(x) x(:)', durations, 'UniformOutput', false);
toplot = [toplot{:}];
histogram(toplot);hold on
text(0.7*max(xlim),0.8*max(ylim),('Median: '+string(median(toplot,'omitnan'))))
title('Distribution of timing between trials for PD with DBS data')
%saveas(gcf,'Timing with median.png')

disp('Total median:   '+string(median(toplot,'omitnan')))
disp('Total trials:   '+string(sum(cellfun(@(x) x,n_trials))))
diary off



