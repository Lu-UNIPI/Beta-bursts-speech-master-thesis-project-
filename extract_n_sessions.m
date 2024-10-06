clear all;close all;clc;
%% Get n session per DBS and ECoG, with and without LFPs
% This script is for data exploration. It calculates n patients with DBS in STN and VIM,
% and n sessions per each cases 

PD=[];
PD.ns_ecog_tot=0;
PD.ns_DBS=0;
PD.ns_ecog_DBS=0;
PD.n_pat_tot=0;
PD.n_pat_DBS=0;
ET=[];
ET.ns_ecog_tot=0;
ET.ns_DBS=0;
ET.ns_ecog_DBS=0;
ET.n_pat_tot=0;
ET.n_pat_DBS=0;
PATH_DATA='Z:\DBS';

n_sub_PD=3001:3032;
n_sub_PD=arrayfun(@(x) sprintf('%04d', x), n_sub_PD, 'UniformOutput', false);
PD.ind_ecogwithDBS=false(1,numel(n_sub_PD));
n_sub_ET=[4057,4058,4060,4062,4066,4067,4069:4078,4080,4083,4084,4086,4087];
n_sub=[n_sub_PD,num2cell(n_sub_ET)];
for i=1:numel(n_sub)
    % open a subject dir
    SUBJECT=strcat('DBS',string(n_sub(i)));
    PATH_SUBJECT_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    
    % retrieve n session annot and type info
    if ~exist(PATH_SUBJECT_ANNOT)
        continue
    end
    T=bml_annot_read(strcat(PATH_SUBJECT_ANNOT, filesep, strcat(SUBJECT, '_session.txt')));
    
    %calculate n_sessions
    if i<=length(n_sub_PD)
        PD.n_pat_tot=PD.n_pat_tot+1;
        PD.ns_ecog_tot=PD.ns_ecog_tot+height(T);
        ns_LMA=cellfun(@(x) sum(ismember(string(T.type),x)),"LMA");
        ns_Lead=cellfun(@(x) sum(ismember(string(T.type),x)),"Lead");
        if (ns_LMA+ns_Lead)>0
            PD.ind_ecogwithDBS(i)=true;
            PD.n_pat_DBS=PD.n_pat_DBS+1;
            PD.ns_ecog_DBS=PD.ns_ecog_DBS+height(T);
            PD.ns_DBS=PD.ns_DBS+ns_LMA+ns_Lead;
        end
    else
        ET.n_pat_tot=ET.n_pat_tot+1;
        ET.ns_ecog_tot=ET.ns_ecog_tot+height(T);
        % ns_LMA=cellfun(@(x) sum(ismember(string(T.type),x)),"LMA"); what is LMA??
        ns_Lead=cellfun(@(x) sum(ismember(string(T.type),x)),"Lead");
        ns_LEAD=cellfun(@(x) sum(ismember(string(T.type),x)),"LEAD");
        if (ns_LMA+ns_Lead+ns_LEAD)>0
            ET.n_pat_DBS=ET.n_pat_DBS+1;
            ET.ns_ecog_DBS=ET.ns_ecog_DBS+height(T);
            ET.ns_DBS=ET.ns_DBS+ns_LMA+ns_Lead+ns_LEAD;
        end
    end
end
%save("PD_data.mat","PD")
%save("ET_data.mat","ET")

%% Visual output
y=[PD.ns_DBS,PD.ns_ecog_DBS,PD.ns_ecog_tot; ET.ns_DBS,ET.ns_ecog_DBS,ET.ns_ecog_tot];
p=bar(y);
legend("DBS", "ECoG with DBS", "ECoG all")
xticks(1:2)%length(y))
x_1=string(PD.n_pat_tot)+"("+string(PD.n_pat_DBS)+") PD patients"+"\newline"+"Electordes in STN";
x_2=string(ET.n_pat_tot)+"("+string(ET.n_pat_DBS)+") ET patients and one Dystonia"+"\newline"+"Electordes in VIM";
xticklabels({x_1,x_2})
title('N sessions per PD and ET data')
for i=1:3
    xtips = p(i).XEndPoints;
    ytips = p(i).YEndPoints;
    labels = string(p(i).YData);
    text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
end
set(gcf,'Position',[300 200 600 500])
%saveas(gcf,'N_sessions.png')

