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
PATH_MNI_BRAIN='C:\...cortex/CortexHiRes.mat';
ld=load(PATH_MNI_BRAIN, 'Vertices', 'Faces');
cortex=reducepatch(ld.Faces, ld.Vertices, 0.3);
faces_cortex = cortex.faces;
vertices_cortex = cortex.vertices;

% data for stn
load("C:\...\atlas_index.mat");
target={'STN','GPi'};
faceColors={[1 0.9 0.7],[0.7 0.9 0.7]};
faceAlphas=[0.4 0.4];
gray=[0.8 0.8 0.8];
atlas_names = cellfun(@(x) x.name, atlases.roi, 'UniformOutput', false);
[~,atlas_idx] = ismember(target, atlas_names);
assert(any(~isempty(atlas_idx)),"No target in this atlas. Change the Atlas!");
faces_bg = arrayfun(@(x) atlases.roi{x,2}.fv.faces,atlas_idx,'UniformOutput',false);
vertices_bg = arrayfun(@(x) atlases.roi{x,2}.fv.vertices,atlas_idx,'UniformOutput',false);

colors=lin_bwr_colormap();


%% analisi
measures={'frequency','probability','duration','volt_amp','time_rdsym','time_ptsym'};%'power'
subdivisions={'Volume','Accuracy','Reaction time','Speech rate','Jitter mean','Jitter CV','Shimmer mean','Shimmer CV'};
for div=1:numel(subdivisions)
subdivision=subdivisions{div};
for m=1:numel(measures)
    fig1=figure('units','normalized','outerposition',[0.03 0.03 0.95 0.95]);
    set(gcf,'Visible','off')

    % -----------------PUT HERE WHAT TO PLOT!!-------------------------
    toEvaluate='fstat_all';
    pv_toEvaluate='pv_all';
    measure=measures{m};
    modality='positive-negative';
    % -----------------------------------------------------------------

    sgtitle(strcat(upper(strrep(measures{m},'_',' ')),"  ",subdivision,"  ",strrep(toEvaluate,'_',' ')))

 
    %% get data 
    ii=1:numel(SUBJECTS);
    alldata=cell(size(ii));

    windows={'baseline','stimulus','prespeech','speech','rebound'}; 
    limits=struct();
    for i=ii
        %open a subject dir
        if i==6;continue;end
        if i==13;continue;end
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
        
        % get BURSTS DATA
        tab=readtable(strcat('annot/...,SUBJECT,'_anova2_',strrep(measure,'_',''),'.txt'));
        alldata{i}=tab;
    end
          
    disp('-----------------------------------')
%     ii=1:numel(SUBJECTS);
%     alldata=cell(size(ii));
% 
%     windows={'baseline','stimulus','prespeech','speech','rebound'}; 
%     limits=struct();
%     for i=ii
%         %open a subject dir
%         if i==6;continue;end
%         if i==14;continue;end
%         SUBJECT=strcat('DBS',string(SUBJECTS(i)));
%         disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
%         
%         % get BURSTS DATA
%         tab=readtable(strcat('annot/...,'_anova2_',strrep(measure,'_',''),'.txt'));
%         alldata{i}=tab(strcmp(tab.trial_division,subdivision),:);
%         tab_ecog=tab(startsWith(tab.label,'ecog'),:);
%         tab_dbs=tab(startsWith(tab.label,'dbs_L'),:);        
%         
%         burst_feature=tab_ecog.(toEvaluate); % array
%         limits.min_ecog(i)=min(burst_feature);%prctile(burst_feature,10,'all');
%         limits.max_ecog(i)=max(burst_feature);%prctile(burst_feature,90,'all');
%         
%         burst_feature=tab_dbs.(toEvaluate); % array
%         limits.min_dbs(i)=min(burst_feature);%prctile(burst_feature,10,'all');
%         limits.max_dbs(i)=max(burst_feature);%prctile(burst_feature,90,'all');
%     
%     end
%     
%     limits.min_ecog = min(limits.min_ecog,[],'all');
%     limits.max_ecog = max(limits.max_ecog,[],'all');
%     limits.min_dbs = min(limits.min_dbs,[],'all');
%     limits.max_dbs = max(limits.max_dbs,[],'all');
%     
%     disp('-----------------------------------')
%     
            
    %% ecog
    subplot(1,2,1)
    cfg = [];
    cfg.h_ax = gca;
    cfg.surface_facealpha = 0.4;
    %cfg.surface_facecolor = [0.7 0.8 1];
    cfg.view = [-100,15];
    bml_plot3d_surface(cfg, vertices_cortex,faces_cortex);hold on
    
    
    for i=ii
        if i==6;continue;end
        if i==13;continue;end
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
        % paths 
        PATH_ANNOT=strcat(PATH_DATA, ...\annot');
        electrode=bml_annot_read(strcat(PATH_ANNOT,...,'_electrode'));
        cfg=[];
        cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                    % threshold_subj or threshold_weight
        electrode=bml_getEcogArea(cfg,electrode);
    
        % load data
        tab=alldata{i};
        tab_ecog=tab(startsWith(tab.label,'ecog'),:);
        channels_ecog=unique(tab_ecog.label); 
        for e=1:numel(channels_ecog)
            if strcmp(electrode.HCPMMP1_area{e},'out of atlas') || strcmp(electrode.HCPMMP1_area{e},'label not present');continue;end
            
            value = tab_ecog.(toEvaluate)( strcmp(tab_ecog.label,channels_ecog(e)) &  strcmp(tab_ecog.trial_division,subdivision) );
            pv = tab_ecog.(pv_toEvaluate)( strcmp(tab_ecog.label,channels_ecog(e)) &  strcmp(tab_ecog.trial_division,subdivision) );
            if pv<0.0000001;pv=0.0000001;end
            
            if i==14;idx=find(strcmp(electrode.electrode(1:130,:),channels_ecog{e}));
                else;idx=find(strcmp(electrode.electrode,channels_ecog{e}));
            end
            point = table();
            point.x = electrode.mni_nonlinear_x(idx);
            point.y = electrode.mni_nonlinear_y(idx);
            point.z = electrode.mni_nonlinear_z(idx);

            if pv<0.05
                cfg=[];
                cfg.limits='standard';
                cfg.modality=modality; % or 'positive'
                cfg.location='ecog';
                value_color=value_to_color(cfg,value);   
                point.color = rgb2hex(colors(round(value_color),:));
                
                cfg=[];
                cfg.limits='standard';
                cfg.modality=modality;
                cfg.location='ecog';
                cfg.estremes=[0.2,3];
                radius=value_to_size(cfg,value);          
            else
                continue
                point.color=gray;
                radius=0.5;
            end
            point.radius = radius;
            point.name = channels_ecog(e);
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
    %title(toEvaluate)
        
        
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
        if i==6;continue;end
        if i==13;continue;end
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
        % paths 
        PATH_ANNOT=strcat(PATH_DATA, ..., 'Preprocessed data\Sync\annot');
        electrode=bml_annot_read(strcat(PATH_ANNOT,...,'_electrode'));
        cfg=[];
        cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                    % threshold_subj or threshold_weight
        electrode=bml_getEcogArea(cfg,electrode);
    
        % load data
        tab=alldata{i};
        tab_dbs=tab(startsWith(tab.label,'dbs_L'),:); 

        % check if there is the need of differentated coords
        if height(tab_dbs)>3
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


        channels_dbs=unique(tab_dbs.label);
        for e=1:numel(channels_dbs)
            value = tab_dbs.(toEvaluate)( strcmp(tab_dbs.label,channels_dbs(e)) &  strcmp(tab_dbs.trial_division,subdivision) );
            pv = tab_dbs.(pv_toEvaluate)( strcmp(tab_dbs.label,channels_dbs(e)) &  strcmp(tab_dbs.trial_division,subdivision) );
            if pv<0.0000001;pv=0.0000001;end
    
            chan=channels_dbs{e};
            minuend=char(chan(1:strfind(chan,'-')-1));
            subtrahend=char(chan(strfind(chan,'-')+1:end));
    
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
           
            if pv<0.05
                cfg=[];
                cfg.limits='standard';
                cfg.location='dbs';
                value_color=value_to_color(cfg,value);   
                point.color = rgb2hex(colors(round(value_color),:));
                
                cfg=[];
                cfg.limits='standard';
                cfg.location='dbs';
                cfg.estremes=[0.1,0.55];
                radius=value_to_size(cfg,value);
            else
                continue
                point.color = gray;
                radius=0.1;
            end
            
            point.radius = radius;
            point.name = channels_dbs(e);
      
            cfg.h_ax = gca;
            cfg.annotate = false;
            bml_plot3d_points(cfg, point);
        end
    colormap(colors);
    c=colorbar;
    % ticks and ticklabels
    if strcmp(modality,'positive');min_value=0;max_value=5;
    else;min_value=-5;max_value=5;
    end
%     min_value = limits.min_dbs;max_value = limits.max_dbs;
    c.Ticks = [0,1]; 
    c.TickLabels = {num2str(min_value), num2str(max_value)};  
    %title(toEvaluate)
    end

    saveas(fig1,strcat('images/bursts CTAR/anova2/',measures{m}," ",subdivision," ",toEvaluate,'.png'));
    saveas(fig1,strcat('images/bursts CTAR/anova2/',measures{m}," ",subdivision," ",toEvaluate,'.fig'));
    close all

end
end

disp(strcat('Finish with: ',toEvaluate))
