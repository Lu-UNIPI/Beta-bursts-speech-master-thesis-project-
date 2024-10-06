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
target={'STN'};
%faceColors={'#F7B05B','#B2E4DB'};
%faceAlphas=[0.6 0.6];
atlas_names = cellfun(@(x) x.name, atlases.roi, 'UniformOutput', false);
[~,atlas_idx] = ismember(target, atlas_names);
assert(any(~isempty(atlas_idx)),"No target in this atlas. Change the Atlas!");
faces_stn = cell2mat( arrayfun(@(x) atlases.roi{x,2}.fv.faces,atlas_idx,'UniformOutput',false));
vertices_stn = cell2mat (arrayfun(@(x) atlases.roi{x,2}.fv.vertices,atlas_idx,'UniformOutput',false));

colors=lin_bwr_colormap();

%% get data
ii=1:length(SUBJECTS);
alldata=cell(size(ii));
windows={'base_beta','stim_beta','prespeech_beta','speech_beta','rebound_beta'};
tab=table();
limits=struct();
for i=ii
    if i==6;continue;end
    %open a subject dir
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    % --------------- PUT HERE WHAT TO PLOT -----------------------------
    tab_pos=readtable(strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_volume_pos.txt')); 
    tab_neg=readtable(strcat('annot/conditions CTAR/',SUBJECT,'_speech_normdb_volume_neg.txt'));
    variables=windows;variables(6)={'id'};variables(7)={'label'};
    tab=tab_pos(:,variables);
    for w=1:numel(windows)
        window=string(windows(w));
        tab.(window)=(tab_pos.(window)-tab_neg.(window))./(tab_pos.(window)+tab_neg.(window));
    end
    alldata{i}=tab;
    for w=1:numel(windows)
        window=string(windows(w));
        disp(strcat('Getting limits :'," ",window))

        pow_ecog=tab.(window)(startsWith(tab.label,'ecog')); % array
        if ~isempty(pow_ecog( pow_ecog<0))
            limits.min_neg_ecog(i,w)=prctile(pow_ecog( pow_ecog<0 ),10,'all');%min( pow_ecog( pow_ecog<0 ) );
            limits.max_neg_ecog(i,w)=prctile(pow_ecog( pow_ecog<0 ),90,'all');%max( pow_ecog( pow_ecog<0 ) );
        end
        if  ~isempty(pow_ecog( pow_ecog>0))
            limits.min_pos_ecog(i,w)=prctile(pow_ecog( pow_ecog>0 ),10,'all');%min( pow_ecog( pow_ecog>0 ) );
            limits.max_pos_ecog(i,w)=prctile(pow_ecog( pow_ecog>0 ),90,'all');%max( pow_ecog( pow_ecog>0 ) );
        end

        pow_dbs=tab.(window)(startsWith(tab.label,'dbs_L'),:); % array
        if ~isempty(pow_dbs( pow_dbs<0))
            limits.min_neg_dbs(i,w)=prctile(pow_dbs( pow_dbs<0 ),10,'all');%min( pow_dbs( pow_dbs<0 ) );
            limits.max_neg_dbs(i,w)=prctile(pow_dbs( pow_dbs<0 ),90,'all');%max( pow_dbs( pow_dbs<0 ) );
        end
        if  ~isempty(pow_dbs( pow_dbs>0))
            limits.min_pos_dbs(i,w)=prctile(pow_dbs( pow_dbs>0 ),10,'all');%min( pow_dbs( pow_dbs>0 ) );
            limits.max_pos_dbs(i,w)=prctile(pow_dbs( pow_dbs>0 ),90,'all');%max( pow_dbs( pow_dbs>0 ) );prctile(pow_dbs( pow_dbs>0 ),90,'all');
        end

    end
end

limits.min_neg_ecog = min(limits.min_neg_ecog,[],'all');
limits.max_neg_ecog = max(limits.max_neg_ecog,[],'all');
limits.min_pos_ecog = min(limits.min_pos_ecog,[],'all');
limits.max_pos_ecog = max(limits.max_pos_ecog,[],'all');
limits.min_neg_dbs = min(limits.min_neg_dbs,[],'all');
limits.max_neg_dbs = max(limits.max_neg_dbs,[],'all');
limits.min_pos_dbs = min(limits.min_pos_dbs,[],'all');
limits.max_pos_dbs = max(limits.max_pos_dbs,[],'all');

disp('-----------------------------------')

%% analisi
figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
% set(gcf,'Visible','off') % 2 views
% t=tiledlayout(2,2);
% title(t,tit_figure)

subplot(2,numel(windows)+1,numel(windows)+1)
axis off; colormap(colors);
c=colorbar;
c.Location = 'west'; 
c.Label.String = 'diff power normalized (dB)'; 
% ticks and ticklabels
min_value = limits.min_neg_ecog;max_value = limits.max_pos_ecog;
c.Ticks = [0,1]; 
c.TickLabels = {num2str(min_value), num2str(max_value)}; 

subplot(2,numel(windows)+1,2*(numel(windows)+1))
axis off; colormap(colors);
c=colorbar;
c.Location = 'west'; 
c.Label.String = 'diff power normalized (dB)'; 
% ticks and ticklabels
min_value = limits.min_neg_dbs;max_value = limits.max_pos_dbs;
c.Ticks = [0,1]; 
c.TickLabels = {num2str(min_value), num2str(max_value)}; 



for w=1:numel(windows)
window=string(windows(w));
disp(strcat('ANALYSIS :'," ",window))  

% nexttile(j)
% cfg = [];
% cfg.h_ax = gca;
% cfg.surface_facealpha = 0.5;
% cfg.view = [-100,15];
% bml_plot3d_surface(cfg, vertices,faces);hold on
% title(strrep(window,'_',' '))


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
    pow_ecog=tab.(window)(startsWith(tab.label,'ecog')); % array
    pow_dbs=tab.(window)(startsWith(tab.label,'dbs_L'),:); % array
    channels_ecog=tab.label(startsWith(tab.label,'ecog'));
    channels_dbs=tab.label(startsWith(tab.label,'dbs_L'));


    %% ecog
    subplot(2,numel(windows)+1,w)
    cfg = [];
    cfg.h_ax = gca;
    cfg.surface_facealpha = 0.5;
    cfg.view = [-100,15];
    bml_plot3d_surface(cfg, vertices_cortex,faces_cortex);hold on

    title(strrep(window,'_',' '))

    min_neg=limits.min_neg_ecog; max_neg=limits.max_neg_ecog;
    min_pos=limits.min_pos_ecog; max_pos=limits.max_pos_ecog;

    for e=1:numel(channels_ecog)
        if strcmp(electrode.HCPMMP1_area{e},'out of atlas') || strcmp(electrode.HCPMMP1_area{e},'label not present');continue;end
        idx=find(strcmp(electrode.electrode,channels_ecog{e}));
        point = table();
        point.x = electrode.mni_nonlinear_x(idx);
        point.y = electrode.mni_nonlinear_y(idx);
        point.z = electrode.mni_nonlinear_z(idx);
        point.radius = 0.5*2;
%         pow_window=tab.(window);
%         if pow_window(e)<0
%             value_to_color = (pow_window(e) - min(pow_window(pow_window<0))) / (max(pow_window(pow_window<0)) - min(pow_window(pow_window<0))) *128;
%         elseif pow_window(e)>0
%             if pow_window(e)==min(pow_window(pow_window>0)) && pow_window(e)==max(pow_window(pow_window>0));value_to_color=128+128/2;
%             else;value_to_color = 128 + ( pow_window(e) - min(pow_window(pow_window>0)) ) / (max(pow_window(pow_window>0)) - min(pow_window(pow_window>0))) *128;
%             end
%         end
        
        if pow_ecog(e)<0
            if pow_ecog(e)==min_neg && pow_ecog(e)==max_neg;value_to_color=128/2; % if there is only a value over zero
            elseif pow_ecog(e)<min_neg;value_to_color=0; % if the value is an outlier
            elseif pow_ecog(e)>max_neg;value_to_color=128;
            else;value_to_color = (pow_ecog(e) - min_neg) / (max_neg - min_neg) *128; 
            end
        elseif pow_ecog(e)>0
            if pow_ecog(e)==min_pos && pow_ecog(e)==max_pos;value_to_color=128+128/2;  % if there is only a value below zero
            elseif pow_ecog(e)>max_pos;value_to_color=256; % if the value is an outlier
            elseif pow_ecog(e)<min_pos;value_to_color=128;
            else;value_to_color = 128 + ( pow_ecog(e) - min_pos ) / (max_pos - min_pos) *128;
            end
        else 
            value_to_color=128;
        end
        if round(value_to_color)==0;value_to_color=1;end
        point.color = rgb2hex(colors(round(value_to_color),:));
        point.name = channels_ecog(e);

        cfg.h_ax = gca;
        cfg.annotate = false;
        bml_plot3d_points(cfg, point);
    end

    %% dbs
    subplot(2,numel(windows)+1,numel(windows)+1+w)
    cfg = [];
    cfg.h_ax = gca;
    cfg.surface_facealpha = 0.6;
    cfg.view = [-100,15];
    bml_plot3d_surface(cfg,vertices_stn,faces_stn);hold on

    min_neg=limits.min_neg_dbs; max_neg=limits.max_neg_dbs;
    min_pos=limits.min_pos_dbs; max_pos=limits.max_pos_dbs;

    for e=1:numel(channels_dbs)
        chan=channels_dbs{e};
        minuend=char(chan(1:strfind(chan,'-')-1));
        subtrahend=char(chan(strfind(chan,'-')+1:end));

        point = table();
        point.x = ( electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_x(find(strcmp(electrode.electrode,subtrahend))) ) /2;
        point.y = ( electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_y(find(strcmp(electrode.electrode,subtrahend))) ) /2;
        point.z = ( electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,minuend))) + electrode.mni_nonlinear_z(find(strcmp(electrode.electrode,subtrahend))) ) /2;
        point.radius = 0.2;

        if pow_dbs(e)<0
            if pow_dbs(e)==min_neg && pow_dbs(e)==max_neg;value_to_color=128/2;
            elseif pow_dbs(e)<min_neg;value_to_color=0; % if the value is an outlier
            elseif pow_dbs(e)>max_neg;value_to_color=128;
            else;value_to_color = (pow_dbs(e) - min_neg) / (max_neg - min_neg) *128;
            end
        elseif pow_dbs(e)>0
            if pow_dbs(e)==min_pos && pow_dbs(e)==max_pos;value_to_color=128+128/2;
            elseif pow_dbs(e)>max_pos;value_to_color=256; % if the value is an outlier
            elseif pow_dbs(e)<min_pos;value_to_color=128;
            else;value_to_color = 128 + ( pow_dbs(e) - min_pos ) / (max_pos - min_pos) *128;
            end
        else 
            value_to_color=128;
        end
        if round(value_to_color)==0;value_to_color=1;end
        point.color = rgb2hex(colors(round(value_to_color),:));
        point.name = channels_dbs(e);

        cfg.h_ax = gca;
        cfg.annotate = false;
        bml_plot3d_points(cfg, point);
    end

end
end

disp('---saving')
saveas(fig1,strcat('images/conditions CTAR/cortex/',SUBJECT,'_cortex_main.png'))
close all

% electrode = table();
% electrode.x = Locmap_error_contrast.X;
% electrode.y = Locmap_error_contrast.Y;
% electrode.z = Locmap_error_contrast.Z;
% electrode.radius = Locmap_error_contrast.Size*2;
% electrode.color = Locmap_error_contrast.Color;%
% % repmat(hex2rgb('#B7B6C1'),height(electrode),1);
%  
%  
% figure; set(gcf,'Visible','off') % 2 views
% tiledlayout(1,2)
% cfg = [];
% cfg.h_ax = gca;
% cfg.view = [-20,15];
% %title({[SUBJECT], ' native-space ECoG reconstruction'})
% hold on
%  
%     %cfg.surface_facecolor = hex2rgb(faceColors{target_i});
%     %cfg.surface_facealpha = faceAlphas(target_i);
%  
% bml_plot3d_surface(cfg, vertices,faces);
%  
%  
% hold on
%  
% cfg.h_ax = gca;
% cfg.annotate = false;
% bml_plot3d_points(cfg, point);
%  
% clim([-3 3])
% colormap(redwhitegreen)
% colorbar