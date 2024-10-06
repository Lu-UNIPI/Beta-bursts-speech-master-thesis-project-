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
faceColors={[1 0.9 0.7],[0.7 0.9 0.7]};
faceAlphas=[0.4 0.4];
gray=[0.8 0.8 0.8];
atlas_names = cellfun(@(x) x.name, atlases.roi, 'UniformOutput', false);
[~,atlas_idx] = ismember(target, atlas_names);
assert(any(~isempty(atlas_idx)),"No target in this atlas. Change the Atlas!");
faces_bg = arrayfun(@(x) atlases.roi{x,2}.fv.faces,atlas_idx,'UniformOutput',false);
vertices_bg = arrayfun(@(x) atlases.roi{x,2}.fv.vertices,atlas_idx,'UniformOutput',false);


%% analisi
measures={'frequency','probability','duration','volt_amp','time_rdsym','time_ptsym'}; %power


for m=1:numel(measures)
    fig1=figure('units','normalized','outerposition',[0.03 0.03 0.95 0.95]);
    set(gcf,'Visible','off')

    % -----------------PUT HERE WHAT TO PLOT!!-------------------------
    toEvaluate='fstat_all';
    pv_toEvaluate='pv_all';
    measure=measures{m};
    modality='positive';
    % -----------------------------------------------------------------
    
    if strcmp(modality,'positive')
        colors = lin_bwr_colormap();
    elseif strcmp(modality,'positive-negative')
        colors = lin_wr_colormap();
    end
    sgtitle(strcat(upper(strrep(measures{m},'_',' ')),"  ",strrep(toEvaluate,'_',' ')))

    %% get data 
    ii=1:numel(SUBJECTS);
    alldata=cell(size(ii));

    windows={'baseline','stimulus','prespeech','speech','rebound'}; 
    limits=struct();
    for i=ii
        %open a subject dir
        if i==6;continue;end
        SUBJECT=strcat('DBS',string(SUBJECTS(i)));
        disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
        
        % get BURSTS DATA
        tab=readtable(strcat('annot/general CTAR/anova1/',SUBJECT,'_anova1_',measure,'.txt'));
        alldata{i}=tab;
    end
          
    disp('-----------------------------------')
%     %% get data 
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
%         tab=readtable(strcat('annot/general CTAR/anova1/',SUBJECT,'_anova1_',strrep(measure,'_',''),'.txt'));
%         alldata{i}=tab;
%         tab_ecog=tab(startsWith(tab.label,'ecog'),:);
%         tab_dbs=tab(startsWith(tab.label,'dbs_L'),:);        
%         
%         burst_feature=tab_ecog.(toEvaluate); % array
%         pv=tab_ecog.(pv_toEvaluate);
%         burst_feature=burst_feature(pv<0.05);
%         if ~isempty(burst_feature(burst_feature<0))
%             limits.min_neg_ecog(i)=min(burst_feature ( burst_feature<0 ));%prctile(burst_feature,10,'all');
%             limits.max_neg_ecog(i)=max(burst_feature ( burst_feature<0 ));%prctile(burst_feature,90,'all');
%         else
%             limits.min_neg_ecog(i)=nan;
%             limits.max_neg_ecog(i)=nan;
%         end
%         if ~isempty(burst_feature(burst_feature>0))
%             limits.min_pos_ecog(i)=min(burst_feature ( burst_feature>0 ));%prctile(burst_feature,10,'all');
%             limits.max_pos_ecog(i)=max(burst_feature ( burst_feature>0 ));%prctile(burst_feature,90,'all');
%         else
%             limits.min_pos_ecog(i)=nan;
%             limits.max_pos_ecog(i)=nan;
%         end
% 
%         burst_feature=tab_dbs.(toEvaluate); % array
%         pv=tab_dbs.(pv_toEvaluate);
%         burst_feature=burst_feature(pv<0.05);
%         if ~isempty(burst_feature(burst_feature<0))
%             %if pv<0.0
%             limits.min_neg_dbs(i)=min(burst_feature ( burst_feature<0 ));%prctile(burst_feature,10,'all');
%             limits.max_neg_dbs(i)=max(burst_feature ( burst_feature<0 ));%prctile(burst_feature,90,'all');
%         else
%             limits.min_neg_dbs(i)=nan;
%             limits.max_neg_dbs(i)=nan;
%         end
%         if ~isempty(burst_feature(burst_feature>0))
%             limits.min_pos_dbs(i)=min(burst_feature ( burst_feature>0 ));%prctile(burst_feature,10,'all');
%             limits.max_pos_dbs(i)=max(burst_feature ( burst_feature>0 ));%prctile(burst_feature,90,'all');
%         else
%             limits.min_pos_dbs(i)=nan;
%             limits.max_pos_dbs(i)=nan;
%         end
%     
%     end
%     
%     limits.min_neg_ecog = min(limits.min_neg_ecog,[],'all');
%     limits.max_neg_ecog = max(limits.max_neg_ecog,[],'all');
%     limits.min_pos_ecog = min(limits.min_pos_ecog,[],'all');
%     limits.max_pos_ecog = max(limits.max_pos_ecog,[],'all');
%     limits.min_neg_dbs = min(limits.min_neg_dbs,[],'all');
%     limits.max_neg_dbs = max(limits.max_neg_dbs,[],'all');
%     limits.min_pos_dbs = min(limits.min_pos_dbs,[],'all');
%     limits.max_pos_dbs = max(limits.max_pos_dbs,[],'all');
%     
%     disp('-----------------------------------')
    
            
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
        tab_ecog=tab(startsWith(tab.label,'ecog'),:);
        %toprint= tab_ecog.(toEvaluate)( tab_ecog.(pv_toEvaluate)<0.05 )
        channels_ecog=unique(tab_ecog.label); 

        for e=1:numel(channels_ecog)
            if strcmp(electrode.HCPMMP1_area{e},'out of atlas') || strcmp(electrode.HCPMMP1_area{e},'label not present');continue;end
            
            value = tab_ecog.(toEvaluate)( strcmp(tab_ecog.label,channels_ecog(e)) );
            pv = tab_ecog.(pv_toEvaluate)( strcmp(tab_ecog.label,channels_ecog(e)) );
            
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

                figure;
                patch([0 1 1 0], [0 0 1 1], colors(round(value_color),:));
                axis off;
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
            value = tab_dbs.(toEvaluate)( strcmp(tab_dbs.label,channels_dbs(e)),: );
            pv = tab_dbs.(pv_toEvaluate)( strcmp(tab_dbs.label,channels_dbs(e)) );
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

    saveas(fig1,strcat('images/bursts CTAR/anova1/',measures{m}," ",toEvaluate,'.png'));
    saveas(fig1,strcat('images/bursts CTAR/anova1/',measures{m}," ",toEvaluate,'.fig'));
    close all

end

disp(strcat('Finish with: ',toEvaluate))