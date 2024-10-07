close all;clear all;clc;
%%
PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

PATH_MNI_BRAIN='C:\Program Files\LeadDBS_Classic\leaddbs\templates\space/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexHiRes.mat';
ld=load(PATH_MNI_BRAIN, 'Vertices', 'Faces');
cortex=reducepatch(ld.Faces, ld.Vertices, 0.3);
faces_cortex =cortex.faces;
vertices_cortex =cortex.vertices;

% data for stn
load("C:\Program Files\LeadDBS_Classic\leaddbs\templates\space\MNI_ICBM_2009b_NLIN_ASYM\atlases\DISTAL Minimal (Ewert 2017)\atlas_index.mat");
target={'STN','GPi'};
faceColors={[1 0.9 0.7],[0.7 0.9 0.7]};
faceAlphas=[0.4 0.4];
gray=[0.8 0.8 0.8];
atlas_names = cellfun(@(x) x.name, atlases.roi, 'UniformOutput', false);
[~,atlas_idx] = ismember(target, atlas_names);
assert(any(~isempty(atlas_idx)),"No target in this atlas. Change the Atlas!");
faces_bg = arrayfun(@(x) atlases.roi{x,2}.fv.faces,atlas_idx,'UniformOutput',false);
vertices_bg = arrayfun(@(x) atlases.roi{x,2}.fv.vertices,atlas_idx,'UniformOutput',false);

colors = lin_bwr_colormap();
modality='positive-negative';
type_electrode='ecog';
%% get data
ii=1:length(SUBJECTS);
alldata=cell(size(ii));
tit_figure='Volume normalized diff - db normalization';
windows={'base_highbeta','stim_highbeta','prespeech_highbeta','speech_highbeta','rebound_highbeta'};
tab=table();
limits=struct();
for i=ii
    if i==6;continue;end
    %open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    % --------------- PUT HERE WHAT TO PLOT -----------------------------
    tab=readtable(strcat('annot/general CTAR/',SUBJECT,'_speech_pow_db.txt')); 
    alldata{i}=tab;
end
disp('-----------------------------------')

%% analisi
for w=1:numel(windows)
    
window=string(windows(w));
disp(strcat('ANALYSIS :'," ",window)) 
if w==1 || w==2 || w==3;continue;end
fig1=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);

if strcmp(type_electrode,'ecog')
    cfg = [];
    cfg.h_ax = gca;
    cfg.surface_facealpha = 0.4;
    cfg.view = [-100,15];
    bml_plot3d_surface(cfg, vertices_cortex,faces_cortex);hold on
else
    for bg=1:2
        cfg = [];
        cfg.h_ax = gca;
        cfg.surface_facealpha = faceAlphas(bg);
        cfg.surface_facecolor = faceColors{bg};
        cfg.view = [80,0];
        bml_plot3d_surface(cfg,vertices_bg{bg},faces_bg{bg});hold on
    end
end

for i=ii
    if i==6;continue;end
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))

    % paths 
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);

    tab=alldata{i}; 
    if strcmp(type_electrode,'ecog')
        pow=tab.(window)(startsWith(tab.label,'ecog')); % array
        channels=tab.label(startsWith(tab.label,'ecog'));
    else
        pow=tab.(window)(startsWith(tab.label,'dbs_L')); % array
        channels=tab.label(startsWith(tab.label,'dbs_L'));
    end
    

    %% ecog
    for e=1:numel(channels)
        value=pow(e);

        if strcmp(type_electrode,'ecog')
            if strcmp(electrode.HCPMMP1_area{e},'out of atlas') || strcmp(electrode.HCPMMP1_area{e},'label not present');continue;end
            idx=find(strcmp(electrode.electrode,channels{e}));
            point = table();
            point.x = electrode.mni_nonlinear_x(idx);
            point.y = electrode.mni_nonlinear_y(idx);
            point.z = electrode.mni_nonlinear_z(idx);
        else
            % check if there is the need of differentated coords
            if length(channels)>3
                coord_original='no';
                radius=0.1;theta= [0 120 240]* pi/180;
                electrodes_coord=electrode(startsWith(electrode.electrode,'dbs_L'),{'electrode','mni_nonlinear_x','mni_nonlinear_y','mni_nonlinear_z'});
                if strcmp(SUBJECT,'DBS3028');electrodes_coord=electrodes_coord(1:end-2,:);end
                L1 = [electrodes_coord.mni_nonlinear_x(1) electrodes_coord.mni_nonlinear_y(1) electrodes_coord.mni_nonlinear_z(1)];
                L2 = [electrodes_coord.mni_nonlinear_x(2) electrodes_coord.mni_nonlinear_y(2) electrodes_coord.mni_nonlinear_z(2)];
                L3 = [electrodes_coord.mni_nonlinear_x(end-1) electrodes_coord.mni_nonlinear_y(end-1) electrodes_coord.mni_nonlinear_z(end-1)];
                L4 = [electrodes_coord.mni_nonlinear_x(end) electrodes_coord.mni_nonlinear_y(end) electrodes_coord.mni_nonlinear_z(end)];
    
                vect= L2 - L1;
                ort_plane= null(vect(:).');
                points=zeros(1,3);
                for p = 1:3
                    points = L2 + radius * cos(theta(p))*ort_plane(:,1)' + sin(theta(p))*ort_plane(:,2)';
                    electrodes_coord.mni_nonlinear_x(1+p)=points(1);
                    electrodes_coord.mni_nonlinear_y(1+p)=points(2);
                    electrodes_coord.mni_nonlinear_z(1+p)=points(3);  
                end
    
                vect= L3 - L4;
                ort_plane= null(vect(:).');
                points=(zeros(1,3));
                for p = 1:3
                    points = L3 + radius * cos(theta(p))*ort_plane(:,1)' + sin(theta(p))*ort_plane(:,2)';
                    electrodes_coord.mni_nonlinear_x(end-p)=points(1);
                    electrodes_coord.mni_nonlinear_y(end-p)=points(2);
                    electrodes_coord.mni_nonlinear_z(end-p)=points(3);  
                end
            else
                coord_original='yes';
            end

            minuend=char(channels{e}(1:strfind(channels{e},'-')-1));
            subtrahend=char(channels{e}(strfind(channels{e},'-')+1:end));
    
            point = table();
            switch coord_original
                case 'yes'
                    point.x = ( electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,subtrahend))) ) /2;
                    point.y = ( electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,subtrahend))) ) /2;
                    point.z = ( electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,subtrahend))) ) /2;
                case 'no'
                    point.x = ( electrodes_coord.mni_nonlinear_x(find(strcmp(electrodes_coord.electrode,minuend))) + electrodes_coord.mni_nonlinear_x(find(strcmp(electrodes_coord.electrode,subtrahend))) ) /2;
                    point.y = ( electrodes_coord.mni_nonlinear_y(find(strcmp(electrodes_coord.electrode,minuend))) + electrodes_coord.mni_nonlinear_y(find(strcmp(electrodes_coord.electrode,subtrahend))) ) /2;
                    point.z = ( electrodes_coord.mni_nonlinear_z(find(strcmp(electrodes_coord.electrode,minuend))) + electrodes_coord.mni_nonlinear_z(find(strcmp(electrodes_coord.electrode,subtrahend))) ) /2;
            end
        end

        cfg=[];
        cfg.limits='standard';
        cfg.modality=modality; % or 'positive'
        cfg.location=type_electrode;
        value_color=value_to_color(cfg,value);   
        point.color = rgb2hex(colors(round(value_color),:));
        
        cfg=[];
        cfg.limits='standard';
        cfg.modality=modality;
        cfg.location=type_electrode;
        if strcmp(type_electrode,'ecog'); cfg.estremes=[0.2,3];else;cfg.estremes=[0.1,0.55];end
        radius=value_to_size(cfg,value);

        point.radius = radius;
        point.name = channels(e);
        cfg.h_ax = gca;
        cfg.annotate = false;
        bml_plot3d_points(cfg, point);
    end
end

colormap(colors);
c=colorbar;
% ticks and ticklabels
if strcmp(modality,'positive');min_value=0;max_value=5;
else;min_value=-5;max_value=5;
end
%min_value = limits.min_ecog;max_value = limits.max_ecog;
c.Ticks = [0,1]; 
c.TickLabels = {num2str(min_value), num2str(max_value)};
title(strrep(window,'_',' '))

end

% saveas(gcf,strcat('images/general CTAR/cortex/cortex_stimulus_whole beta.png'))
% saveas(gcf,strcat('images/general CTAR/cortex/cortex_stimulus_whole beta.fig'))
% %% DBS!
% 
% for w=1:numel(windows)
%     
% window=string(windows(w));
% disp(strcat('ANALYSIS :'," ",window)) 
% if w==1 || w==3;continue;end
% fig1=figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
% cfg = [];
% cfg.h_ax = gca;
% cfg.surface_facealpha = 0.4;
% cfg.view = [-100,15];
% bml_plot3d_surface(cfg, vertices_cortex,faces_cortex);hold on
% 
% for i=ii
%     if i==6;continue;end
%     SUBJECT=strcat('DBS',string(SUBJECTS(i)));
%     disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
% 
% % paths 
% PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
% electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
% cfg=[];
% cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
%                             % threshold_subj or threshold_weight
% electrode=bml_getEcogArea(cfg,electrode);
% 
% tab=alldata{i}; 
% pow_dbs=tab.(window)(startsWith(tab.label,'dbs_L'),:); % array
% channels_dbs=tab.label(startsWith(tab.label,'dbs_L'));
% 
%     %% dbs
%     subplot(2,numel(windows)+1,numel(windows)+1+w)
%     cfg = [];
%     cfg.h_ax = gca;
%     cfg.surface_facealpha = 0.4;
%     cfg.view = [-100,15];
%     bml_plot3d_surface(cfg,vertices_stn,faces_stn);hold on
% 
%     min_neg=limits.min_neg_dbs; max_neg=limits.max_neg_dbs;
%     min_pos=limits.min_pos_dbs; max_pos=limits.max_pos_dbs;
% 
%     for e=1:numel(channels_dbs)
%         chan=channels_dbs{e};
%         minuend=char(chan(1:strfind(chan,'-')-1));
%         subtrahend=char(chan(strfind(chan,'-')+1:end));
% 
%         point = table();
%         point.x = ( electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,subtrahend))) ) /2;
%         point.y = ( electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,subtrahend))) ) /2;
%         point.z = ( electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,subtrahend))) ) /2;
%         point.radius = 0.2;
% 
%         if pow_dbs(e)<0
%             if pow_dbs(e)==min_neg && pow_dbs(e)==max_neg;value_to_color=128/2;
%             elseif pow_dbs(e)<min_neg;value_to_color=0; % if the value is an outlier
%             elseif pow_dbs(e)>max_neg;value_to_color=128;
%             else;value_to_color = (pow_dbs(e) - min_neg) / (max_neg - min_neg) *128;
%             end
%         elseif pow_dbs(e)>0
%             if pow_dbs(e)==min_pos && pow_dbs(e)==max_pos;value_to_color=128+128/2;
%             elseif pow_dbs(e)>max_pos;value_to_color=256; % if the value is an outlier
%             elseif pow_dbs(e)<min_pos;value_to_color=128;
%             else;value_to_color = 128 + ( pow_dbs(e) - min_pos ) / (max_pos - min_pos) *128;
%             end
%         else 
%             value_to_color=128;
%         end
%         if round(value_to_color)==0;value_to_color=1;end
%         point.color = rgb2hex(colors(round(value_to_color),:));
%         point.name = channels_dbs(e);
% 
%         cfg.h_ax = gca;
%         cfg.annotate = false;
%         bml_plot3d_points(cfg, point);
%     end
% 
% 
% end
% 
% end
% 
% disp('---saving')
% saveas(gcf,strcat('images/general CTAR/cortex/',SUBJECT,'_cortex_main.png'))
% saveas(gcf,strcat('images/general CTAR/cortex/',SUBJECT,'_cortex_main.fig'))
% close all
