close all;clear all;clc;
%%
PATH_DATA='Z:\DBS';

DATE=datestr(now,'yyyymmdd');
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

% data for cortex
PATH_MNI_BRAIN='C:\Program Files\LeadDBS_Classic\leaddbs\templates\space/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexHiRes.mat';
ld=load(PATH_MNI_BRAIN, 'Vertices', 'Faces');
cortex=reducepatch(ld.Faces, ld.Vertices, 0.3);
faces_cortex = cortex.faces;
vertices_cortex = cortex.vertices;

% data for stn
load("C:\Program Files\LeadDBS_Classic\leaddbs\templates\space\MNI_ICBM_2009b_NLIN_ASYM\atlases\DISTAL Minimal (Ewert 2017)\atlas_index.mat");
target={'STN','GPi'};
faceColors={[1 0.8 0.2],[0.5 0.8 0.5]};
faceAlphas=[0.4 0.4];
atlas_names = cellfun(@(x) x.name, atlases.roi, 'UniformOutput', false);
[~,atlas_idx] = ismember(target, atlas_names);
assert(any(~isempty(atlas_idx)),"No target in this atlas. Change the Atlas!");
faces_bg = arrayfun(@(x) atlases.roi{x,2}.fv.faces,atlas_idx,'UniformOutput',false);
vertices_bg = arrayfun(@(x) atlases.roi{x,2}.fv.vertices,atlas_idx,'UniformOutput',false);

colors=lin_bwr_colormap();
%colors=parula();

%% get data 
ii=1:7;%:4;
alldata=cell(size(ii));

% -----------------PUT HERE WHAT TO PLOT!!-------------------------
toEvaluate='duration';
% -----------------------------------------------------------------

windows={'baseline','stimulus','prespeech','speech','rebound'}; 
limits=struct();
for i=ii
    %open a subject dir
    if i==6;continue;end
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
    % get BURSTS DATA
    tab=readtable(strcat('annot/general CTAR/',SUBJECT,'_beta_bursts.txt')); tab=tab(tab.duration > 0.1, :);
    alldata{i}=tab;
    tab_ecog=tab(startsWith(tab.chan_id,'ecog'),:);
    tab_dbs=tab(startsWith(tab.chan_id,'dbs_L'),:);
    
    for w=1:numel(windows)
        window=string(windows(w));
        disp(strcat('Getting limits :'," ",window))

        burst_feature=tab_ecog.(toEvaluate)(strcmp(tab_ecog.window,window)); % array
        limits.min_ecog(i,w)=prctile(burst_feature,10,'all');
        limits.max_ecog(i,w)=prctile(burst_feature,90,'all');
        
        burst_feature=tab_dbs.(toEvaluate)(strcmp(tab_dbs.window,window)); % array
        limits.min_dbs(i,w)=prctile(burst_feature,10,'all');
        limits.max_dbs(i,w)=prctile(burst_feature,90,'all');
    end
end

limits.min_ecog = min(limits.min_ecog,[],'all');
limits.max_ecog = max(limits.max_ecog,[],'all');
limits.min_dbs = min(limits.min_dbs,[],'all');
limits.max_dbs = max(limits.max_dbs,[],'all');

disp('-----------------------------------')

%% analisi

tit_figure=toEvaluate;

%figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
%set(gcf,'Visible','off') % 2 views
% t=tiledlayout(2,4);
% title(t,tit_figure)

% subplot(2,numel(windows)+1,numel(windows)+1)
% axis off; colormap(colors);
% c=colorbar;
% c.Location = 'west'; 
% c.Label.String = toEvaluate; 
% % ticks and ticklabels
% min_value = limits.min_ecog;max_value = limits.max_ecog;
% c.Ticks = [0,1]; 
% c.TickLabels = {num2str(min_value), num2str(max_value)}; 
% 
% subplot(2,numel(windows)+1,2*(numel(windows)+1))
% axis off; colormap(colors);
% c=colorbar;
% c.Location = 'west'; 
% c.Label.String = toEvaluate; 
% % ticks and ticklabels
% min_value = limits.min_dbs;max_value = limits.max_dbs;
% c.Ticks = [0,1]; 
% c.TickLabels = {num2str(min_value), num2str(max_value)}; 

for w=3:numel(windows)
    window=string(windows(w));
    disp(strcat('ANALYSIS :'," ",window))  
    
    % nexttile(j)
    % cfg = [];
    % cfg.h_ax = gca;
    % cfg.surface_facealpha = 0.5;
    % cfg.view = [-100,15];
    % bml_plot3d_surface(cfg, vertices,faces);hold on
    % title(strrep(window,'_',' '))
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    sgtitle(strcat(window," ",toEvaluate))
     
    
    %% ecog
    %subplot(2,numel(windows)+1,w)    
    subplot(1,2,1)
    cfg = [];
    cfg.h_ax = gca;
    cfg.surface_facealpha = 0.4;
    %cfg.surface_facecolor = [0.7 0.8 1];
    cfg.view = [-100,15];
    bml_plot3d_surface(cfg, vertices_cortex,faces_cortex);hold on
    
    title(strrep(window,'_',' '))
    for i=ii
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
        % paths 
        PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
        electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
        cfg=[];
        cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                    % threshold_subj or threshold_weight
        electrode=bml_getEcogArea(cfg,electrode);
    
        % load data
        tab=alldata{i};
        toPlot=tab(strcmp(tab.window,window),:);
        toPlot_ecog=toPlot(startsWith(toPlot.chan_id,'ecog'),:);
        channels_ecog=unique(toPlot_ecog.chan_id); 
        for e=1:numel(channels_ecog)
            if strcmp(electrode.HCPMMP1_area{e},'out of atlas') || strcmp(electrode.HCPMMP1_area{e},'label not present');continue;end
            
            value = toPlot_ecog.(toEvaluate)( strcmp(toPlot_ecog.chan_id,channels_ecog(e)),: );
            value= mean (value,'omitnan');
            
            idx=find(strcmp(electrode.electrode,channels_ecog{e}));
            point = table();
            point.x = electrode.mni_nonlinear_x(idx);
            point.y = electrode.mni_nonlinear_y(idx);
            point.z = electrode.mni_nonlinear_z(idx);
            
            if value >= limits.max_ecog
                value_to_color = 256;
            elseif value <= limits.min_ecog
                value_to_color = 1;
            else
                value_to_color= (value - limits.min_ecog) / (limits.max_ecog - limits.min_ecog) * 256;
            end
            if round(value_to_color)==0;value_to_color=1;end
            point.color = rgb2hex(colors(round(value_to_color),:));
            radius= 0.7 + ((value_to_color - 1) / (256 - 1)) * (5 - 0.7);
            point.radius = 0.7;
            point.name = channels_ecog(e);
        
            cfg.h_ax = gca;
            cfg.annotate = false;
            bml_plot3d_points(cfg, point);
        end
    end
    colormap(colors);
    c=colorbar;
    % ticks and ticklabels
    min_value = limits.min_ecog;max_value = limits.max_ecog;
    c.Ticks = [0,1]; 
    c.TickLabels = {num2str(min_value), num2str(max_value)};
    
    
    %% dbs
    %subplot(2,numel(windows)+1,(numel(windows)+1)+w)
    % per stn e per gpi, cicla su facealpha e facecolor
    subplot(1,2,2)
    for bg=1:2
        cfg = [];
        cfg.h_ax = gca;
        cfg.surface_facealpha = faceAlphas(bg);
        cfg.surface_facecolor = faceColors{bg};
        cfg.view = [80,0];
        bml_plot3d_surface(cfg,vertices_bg{bg},faces_bg{bg});hold on
    end
    for i=ii
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
        % paths 
        PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
        electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
        cfg=[];
        cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                    % threshold_subj or threshold_weight
        electrode=bml_getEcogArea(cfg,electrode);
    
        % load data
        tab=alldata{i};
        toPlot=tab(strcmp(tab.window,window),:);
        toPlot_dbs=toPlot(startsWith(toPlot.chan_id,'dbs_L'),:); 
        channels_dbs=unique(toPlot_dbs.chan_id);
        for e=1:numel(channels_dbs)
            value = toPlot_dbs.(toEvaluate)( strcmp(toPlot_dbs.chan_id,channels_dbs(e)),: );
            value=mean(value,'omitnan');
    
            chan=channels_dbs{e};
            minuend=char(chan(1:strfind(chan,'-')-1));
            subtrahend=char(chan(strfind(chan,'-')+1:end));
    
            point = table();
            point.x = ( electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,subtrahend))) ) /2;
            point.y = ( electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,subtrahend))) ) /2;
            point.z = ( electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,subtrahend))) ) /2;
            point.radius = 0.1;
            if value >= limits.max_dbs
                value_to_color = 256;
            elseif value <= limits.min_dbs
                value_to_color = 1;
            else
                value_to_color = (value - limits.min_dbs) / (limits.max_dbs - limits.min_dbs) * 256 ;
            end
            if round(value_to_color)==0;value_to_color=1;end
            point.color = rgb2hex(colors(round(value_to_color),:));
            point.name = channels_dbs(e);
      
            cfg.h_ax = gca;
            cfg.annotate = false;
            bml_plot3d_points(cfg, point);
        end
    colormap(colors);
    c=colorbar;
    % ticks and ticklabels
    min_value = limits.min_dbs;max_value = limits.max_dbs;
    c.Ticks = [0,1]; 
    c.TickLabels = {num2str(min_value), num2str(max_value)};   
    end

end
