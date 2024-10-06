clear all;close all;clc;
%% Get positions of ECoG, with and without LFPs
% This script is for data exploration. It decodes the total coverage of the ecog selected 
% following the HCPMMP1 template, and returns a histogram, bin=roi.

%decide path and thresholds
PATH_DATA='Z:\DBS';
TH=5;
TH_W=0.15;

%create color table
regions=["PMC","POC","IFC","SMC","TPOJ","LTC","DLPFC","IFOC","AAC","IPC"]';
colors={[1,0.7,0],[0,0,0.45],[1,1,0],[0.5,0.9,0.5],[0.3,0,0.5],[0,0.5,0],[0.5,0,0.2],[0,0.8,1],[1,0.1,0.1],[0,0.5,0.9]}';
colors_label=["light orange","dark blue","yellow","light green","dark violet",...
    "dark green","dark red","light blue","red","middle blue"]';
colortab=table(regions,colors,colors_label);


%inizialize struct 
ecog_pos=[];
ecog_pos.tot_sessions=[];
ecog_pos.thp_sessions=[]; %more than th in each patient
ecog_pos.thw_sessions=[]; %with a HCPMMP1_weight_1 > th_w
fns = fieldnames(ecog_pos);
for i=1:length(regions)
    for j=1:length(fns)
        ecog_pos.(fns{j}).(regions(i))=[];
        ecog_pos.(fns{j}).(regions(i)).value=0;
        ecog_pos.(fns{j}).(regions(i)).ind=[];
        ecog_pos.(fns{j}).(regions(i)).ind.name="";
        ecog_pos.(fns{j}).(regions(i)).ind.indexes=[];
    end
end
others=[""];

%subjects in analysis
n_sub_PD=3001:3032;
n_sub_PD=arrayfun(@(x) sprintf('%04d', x), n_sub_PD, 'UniformOutput', false);
n_sub_ET=[4057,4058,4060,4062,4066,4067,4069:4078,4080,4083,4084,4086,4087];
n_sub=[n_sub_PD,num2cell(n_sub_ET)];
for i=1:numel(n_sub)
    % open a subject dir
    SUBJECT=strcat('DBS',string(n_sub(i)));
    PATH_SUBJECT_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');

    % retrieve n eletrod annot and type info
    if ~exist(PATH_SUBJECT_ANNOT)
        continue
    end
    T=bml_annot_read(strcat(PATH_SUBJECT_ANNOT, filesep, strcat(SUBJECT, '_electrode.txt')));
    
    % decoding main probable areas
    % subject analysis

    %initialize containers for ind and values
    PMC=0;POC=0;IFC=0;SMC=0;TPOJ=0;LTC=0;DLPFC=0;IFOC=0;AAC=0;IPC=0;
    PMC_thw=0;POC_thw=0;IFC_thw=0;SMC_thw=0;TPOJ_thw=0;LTC_thw=0;DLPFC_thw=0;IFOC_thw=0;AAC_thw=0;IPC_thw=0;
    PMCi=false(1,height(T));POCi=false(1,height(T));IFCi=false(1,height(T));SMCi=false(1,height(T));TPOJi=false(1,height(T));
    LTCi=false(1,height(T));DLPFCi=false(1,height(T));IFOCi=false(1,height(T));AACi=false(1,height(T));IPCi=false(1,height(T));
    PMCi_thw=false(1,height(T));POCi_thw=false(1,height(T));IFCi_thw=false(1,height(T));SMCi_thw=false(1,height(T));TPOJi_thw=false(1,height(T));
    LTCi_thw=false(1,height(T));DLPFCi_thw=false(1,height(T));IFOCi_thw=false(1,height(T));AACi_thw=false(1,height(T));IPCi_thw=false(1,height(T));
    
    for j=1:height(T)
        label=string(T.HCPMMP1_label_1(j));
        weight=T.HCPMMP1_weight_1(j);
        switch label
        case {'6a','6d','FEF','55b','PEF','6r','6v'}
            PMC=PMC+1;PMCi(j)=true;
            if weight>=TH_W;PMC_thw=PMC_thw+1;PMCi_thw(j)=true;end
        case {'OP1','OP2-3','OP4','PFcm','FOP1','43'}
            POC=POC+1;POCi(j)=true;
            if weight>=TH_W;POC_thw=POC_thw+1;POCi_thw(j)=true;end
        case {'44','45','47l','p47r','IFSa','IFSp','IFJa','IFJp'}
            IFC=IFC+1;IFCi(j)=true;
            if weight>=TH_W;IFC_thw=IFC_thw+1;IFCi_thw(j)=true;end
        case{'1','2','3a','3b','4'}
            SMC=SMC+1;SMCi(j)=true;
            if weight>=TH_W;SMC_thw=SMC_thw+1;SMCi_thw(j)=true;end
        case{'PSL','STV','TPOJ1','TPOJ2','TPOJ3'}
            TPOJ=TPOJ+1;TPOJi(j)=true;
            if weight>=TH_W;TPOJ_thw=TPOJ_thw+1;TPOJi_thw(j)=true;end
        case{'TGd','TE1a','TE2a','TE1m','TE1p','TE2p','PHT','TGV','TF'}
            LTC=LTC+1;LTCi(j)=true;
            if weight>=TH_W;LTC_thw=LTC_thw+1;LTCi_thw(j);end
        case{'9a','9p','9-46d','46','p9-46v','a9-46v','8C','8Av','8Ad','i6-8','s6-8','SFL','8BL'}
            DLPFC=DLPFC+1;DLPFCi(j)=true;
            if weight>=TH_W;DLPFC_thw=DLPFC_thw+1;DLPFCi_thw(j)=true;end
        case{'AVI','MI','AAIC','Pir','PI','FOP2','FOP3','FOP4','FOP5','lg','52','Pol1','Pol2'}
            IFOC=IFOC+1;IFOCi(j)=true;
            if weight>=TH_W;IFOC_thw=IFOC_thw+1;IFOCi_thw(j)=true;end
        case{'STGa','STSda','STSva','STSvp','STSdp','TA2','A4','A5'}
            AAC=AAC+1;AACi(j)=true;
            if weight>=TH_W;AAC_thw=AAC_thw+1;AACi_thw(j)=true;end
        case{'PFop','PF','IP0','IP1','IP2','PFm','PGi','PGs','PFt','PGp'}
            IPC=IPC+1;IPCi(j)=true;
            if weight>=TH_W;IPC_thw=IPC_thw+1;IPCi_thw(j)=true;end
  
        otherwise
            if label~=""
                if ~ismember(label,others);others=[others,label];warning("Unassign label: "+label);end
            end
        end

    end %out of the single subject
    %update 
    all=[PMC,POC,IFC,SMC,TPOJ,LTC,DLPFC,IFOC,AAC,IPC];
    all_thw=[PMC_thw,POC_thw,IFC_thw,SMC_thw,TPOJ_thw,LTC_thw,DLPFC_thw,IFOC_thw,AAC_thw,IPC_thw];
    alli=[PMCi',POCi',IFCi',SMCi',TPOJi',LTCi',DLPFCi',IFOCi',AACi',IPCi'];
    alli_thw=[PMCi_thw',POCi_thw',IFCi_thw',SMCi_thw',TPOJi_thw',LTCi_thw',DLPFCi_thw',IFOCi_thw',AACi_thw',IPCi_thw'];
    for k=1:length(regions)
        %update values
        ecog_pos.tot_sessions.(regions(k)).value= ecog_pos.tot_sessions.(regions(k)).value+all(k);
        ecog_pos.thw_sessions.(regions(k)).value= ecog_pos.thw_sessions.(regions(k)).value+all_thw(k);
        if all(k)>TH;ecog_pos.thp_sessions.(regions(k)).value= ecog_pos.thp_sessions.(regions(k)).value+all(k);end
        %update indexes
        ecog_pos.tot_sessions.(regions(k)).ind.indexes=alli(:,k);
        ecog_pos.thw_sessions.(regions(k)).ind.indexes=+alli_thw(:,k);
        if all(k)>TH;ecog_pos.thp_sessions.(regions(k)).ind.indexes=all(:,k);end
        %update names
        ecog_pos.tot_sessions.(regions(k)).ind.name= string(n_sub(i));
        ecog_pos.thw_sessions.(regions(k)).ind.name= string(n_sub(i));
        if all(k)>TH;ecog_pos.thp_sessions.(regions(k)).ind.name= string(n_sub(i));end
    end
end
others=others(2:end);

%% visual rendering
names=cellstr(regions)';
y=arrayfun(@(x) ecog_pos.tot_sessions.(regions(x)).value,1:numel(names)); %y--> tot
z=arrayfun(@(x) ecog_pos.thp_sessions.(regions(x)).value,1:numel(names)); %z-->thp
w=arrayfun(@(x) ecog_pos.thw_sessions.(regions(x)).value,1:numel(names)); %w-->thw
c=colortab.colors;
for i=1:numel(names) %chose what to plot
    toplot=y(i);
    %toplot=[y(i) w(i)];
    p=bar(i,toplot,'FaceColor',c{i});
    hold on
    for j=1:length(toplot)
        xtips = p(j).XEndPoints;
        ytips = p(j).YEndPoints;
        labels = string(p(j).YData);
        text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
    end
end
title("ECoG coverage")

xticks(1:numel(names))
xticklabels(names)

set(gcf,'Position',[300 100 950 580]);
legend("Premotor cortex","Posterior Opercular Cortex","Inferior Frontal Cortex",...
    "Somatosensory and Motor Cortex","Temporo-Parieto-Occipial Junction","Lateral Temporal Cortex",...
    "DorsoLateral Prefrontal Cortex","Insular and Frontal Opercular Cortex",...
    "Auditory Association Cortex","Inferior Parietal Cortex");
%a="ECoG coverage for w over "+string(th_w)+".png";
%saveas(gcf,a)