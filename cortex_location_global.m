close all;clear all;clc;
%%
% This script performs the following tasks:
% 
% Initializes the workspace and paths.
% Reads and processes data related to brain areas and subject-specific electrode placements.
% Plots 3D brain cortex surfaces and electrode positions.
% Generates a histogram representing the distribution of electrode coverage across different brain areas.

%%
% PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long

tab=readtable("HCPMMP1toAreas.txt");
areas=unique(tab.area);
tab_areas=struct();
for a=1:numel(areas)
    label_area=areas{a};
    tab_areas.(label_area)=0;
end

%% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

% data for cortex
PATH_MNI_BRAIN='C:\Program Files\LeadDBS_Classic\leaddbs\templates\space/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexHiRes.mat';
ld=load(PATH_MNI_BRAIN, 'Vertices', 'Faces');
cortex=reducepatch(ld.Faces, ld.Vertices, 0.3);
faces = cortex.faces;
vertices = cortex.vertices;

figure; set(gcf,'Visible','off') 
tiledlayout(1,2)
cfg = [];
cfg.h_ax = gca;
cfg.surface_facealpha = 0.3;
cfg.view = [-100,15];
bml_plot3d_surface(cfg, vertices,faces);hold on

for i=1:numel(SUBJECTS)
    if i==6;continue;end
    if i==14;continue;end
    %open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    % paths 
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    tab_subject=readtable(strcat("annot\general CTAR\",SUBJECT,'_speech_pow_db.txt'));
    tab_channels=tab_subject.label(startsWith(tab_subject.label,'ecog'));

    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);
    chan_ecog=electrode.electrode(startsWith(electrode.electrode,'ecog'));

    nElec=numel(chan_ecog);    
    for e=1:nElec
        ee=find(strcmp(electrode.electrode,chan_ecog{e}));
        if strcmp([electrode.HCPMMP1_area{ee}],'out of atlas') || strcmp([electrode.HCPMMP1_area{ee}],'label not present');continue;end
        if ~any(strcmp(tab_channels,chan_ecog{e}));continue;end
        
        point = table();
        point.x = electrode.mni_nonlinear_x(ee);
        point.y = electrode.mni_nonlinear_y(ee);
        point.z = electrode.mni_nonlinear_z(ee);
        point.radius = 1.2;
        point.color = electrode.HCPMMP1_color{ee};
        point.name = chan_ecog{e};

        cfg.h_ax = gca;
        cfg.annotate = false;
        bml_plot3d_points(cfg, point);

        % for histogram
        tab_areas.(electrode.HCPMMP1_area{ee})=tab_areas.(electrode.HCPMMP1_area{ee})+1;

    end
end

%% histogram

field=fieldnames(tab_areas);
values=(struct2array(tab_areas));
colors=unique(tab.color);

rgbColors = cellfun(@hex2rgb, colors, 'UniformOutput', false);
rgbColors = vertcat(rgbColors{:});

figure;
b=barh(values);

areas_longname=cell(1,numel(field));
b.FaceColor = 'flat';
for col=1:numel(field)
    color=tab.color(strcmp(tab.area,field(col)));
    color=hex2rgb(color{1});
    b.CData(col,:) = color;
    
    fullname=tab.area_fullname(strcmp(tab.area,field(col)));
    areas_longname{col}=fullname{1};
end

set(gca, 'YTickLabel', areas_longname, 'YTick', 1:numel(areas_longname));
set(gca, 'FontSize', 16);

title('ECoG coverage');
xlabel('# channel');
