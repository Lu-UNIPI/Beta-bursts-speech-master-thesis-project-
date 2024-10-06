close all;clear all;clc;
%%
PATH_DATA='Z:\DBS';
DATE=datestr(now,'yyyymmdd');
format long

tab=readtable("HCPMMP1toAreas.txt");
tab_freq=struct();
areas=unique(tab.area);areas{11}='dbs';
for a=1:numel(areas)
    if a== numel(areas);label_area='dbs';
    else;label_area=areas{a};
    end
    bycycle_freq=strcat(label_area,'_bycycle');
    fooof_freq=strcat(label_area,'_fooof');
    tab_freq.(bycycle_freq)=[];
    tab_freq.(fooof_freq)=[];
end

% choose data to analyze
n_sub_PD_DBS=[3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS=arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS=n_sub_PD_DBS;

%% 
ii=1:numel(SUBJECTS);
j=1;jj=1;
x_var_ecog=[];y_absdiff_ecog=[];y_zscore_ecog=[];
x_var_dbs=[];y_absdiff_dbs=[];y_zscore_dbs=[];
colors=zeros(1,3);
x_f_ecog=[];y_b_ecog=[];
x_f_dbs=[];y_b_dbs=[];
for i=ii
    %open a subject dir
    if i==6;continue;end
    SUBJECT=strcat('DBS',string(SUBJECTS(i)));
    disp(strcat('Now running i= ',string(i),'   aka: ',SUBJECT))
    
    PATH_ANNOT=strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');

    % take sessions with DBS LFP data
    session=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_session'));
    id_session=session.id(strcmpi(session.type, 'LEAD'));

    % upload useful annot tables
    coding=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_coding')); coding=coding(coding.session_id==id_session,:);
    cue=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_cue_precise')); cue=cue(cue.session_id==id_session, :);
    prod_triplet=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_produced_triplet'));prod_triplet=prod_triplet(prod_triplet.session_id==id_session,:);
    electrode=bml_annot_read(strcat(PATH_ANNOT,filesep,SUBJECT,'_electrode'));
    electrode=electrode(:, {'id', 'starts','ends','duration','electrode','connector','port','HCPMMP1_label_1','HCPMMP1_weight_1'});

    cfg=[];
    cfg.decodingtype='basic';   % 'basic', 'bysubject', 'weight'
                                % threshold_subj or threshold_weight
    electrode=bml_getEcogArea(cfg,electrode);

    tab_fooof=readtable(strcat('annot/general CTAR/',SUBJECT,'_fooof_bursts.txt'));
    tab_cycle=readtable(strcat('annot/general CTAR/',SUBJECT,'_cycle_bursts.txt'));

    % filter out bursts <=100ms NB CAMBIA MOLTO SE SI VEDE SOLO SOTTO 100ms!!
    tab_fooof=tab_fooof(tab_fooof.duration>0.1001,:);
    tab_cycle=tab_cycle(tab_cycle.duration>0.1001,:);
    
    chan_selected=unique(tab_cycle.chan_id(~(startsWith(tab_cycle.chan_id,'dbs_R'))));
    for e=1:numel(chan_selected)
        if i==14;label_area=electrode.HCPMMP1_area(strcmp(electrode.electrode(1:130,:),chan_selected{e}));
        else;label_area=electrode.HCPMMP1_area(strcmp(electrode.electrode,chan_selected{e}));
        end
        if isempty(label_area)
            label_area='dbs';
            by_freq=( tab_cycle.frequency(startsWith(tab_cycle.chan_id,chan_selected{e})) )';
            f_freq=tab_fooof.ifreq (startsWith(tab_fooof.chan_id,chan_selected{e}));f_freq=f_freq(1);

            x_var_dbs(jj) = var(by_freq);
            y_absdiff_dbs(jj) = abs(f_freq - mean(by_freq));
            y_zscore_dbs(jj) = (f_freq - mean(by_freq))/std(by_freq);
            x_f_dbs(jj) = f_freq;
            y_b_dbs(jj) = mean(by_freq);
            jj=jj+1;
                
        else 
            if any(strcmp({'out of atlas','outside','label not present'},label_area));continue;end

            color=tab.color(strcmp(tab.area,label_area));
            color=hex2rgb(color{1});
    
            by_freq=( tab_cycle.frequency(startsWith(tab_cycle.chan_id,chan_selected{e})) )';
            f_freq=tab_fooof.ifreq (startsWith(tab_fooof.chan_id,chan_selected{e}));f_freq=f_freq(1);
            
            x_var_ecog(j) = var(by_freq);
            y_absdiff_ecog(j) = abs(f_freq - mean(by_freq));
            y_zscore_ecog(j) = (f_freq - mean(by_freq))/std(by_freq);
            x_f_ecog(j) = f_freq;
            y_b_ecog(j) = mean(by_freq);
            colors(j,:)=color;
            j=j+1;
        end
        bycycle_freq=char(strcat(label_area,'_bycycle'));
        fooof_freq=char(strcat(label_area,'_fooof'));
        tab_freq.(bycycle_freq)=[tab_freq.(bycycle_freq) by_freq];
        tab_freq.(fooof_freq)=[tab_freq.(fooof_freq) f_freq];
    end
end
        
%% fooof freq vs cycles 
fig1=figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle('Comparison fooof and bycycle frequencies per area')
fields=fieldnames(tab_freq);
aa=1;
for a=1:2:numel(fields)
    subplot(2,6,aa)
    label_bycycle=fields{a};
    label_fooof=fields{a+1};
    area_short=char(label_bycycle(1:strfind(label_bycycle,'_')-1));
    if a==numel(fields)-1
        area_long='DBS left';
        color=[0.4,0.6,1];
    else
        area_long=tab.area_fullname(strcmp(tab.area,area_short));
        area_long=area_long{1};
        color=tab.color(strcmp(tab.area,area_short));
        color=hex2rgb(color{1});
    end
    
    histogram(tab_freq.(label_bycycle),40,'FaceColor',color);hold on
    iFreq=mean(tab_freq.(label_fooof));
    xline([iFreq iFreq],'LineWidth',2,'Color',[1,0.7,0.4])
    xrange = [iFreq - 2.5, iFreq + 2.5];
    x_fill = [xrange(1), xrange(1), xrange(2), xrange(2)];
    y_fill = [min(ylim), max(ylim), max(ylim), min(ylim)];
    fill(x_fill, y_fill, [1,0.7,0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    x_lim=xlim;
    y_lim=ylim;
    zscore=(iFreq - mean(tab_freq.(label_bycycle)))/std(tab_freq.(label_bycycle));
    text(x_lim(1) + 0.6*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)), sprintf('Zscr = %.3f \nVar= %.3f',...
            zscore,var(tab_freq.(label_bycycle))),'FontSize',9, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
%     var_index = var(tab_freq.(label_bycycle)) / abs(mean(tab_freq.(label_bycycle))-iFreq);
%     x_lim=xlim;
%     y_lim=ylim;
%     text(x_lim(1) + 0.7*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)), sprintf('Zscore = %.2f \nVar = %.2f  \nVar index = %.2f',...
%             zscore,var(tab_freq.(label_bycycle)),var_index),...
%             'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
    
    title(area_long,'FontSize',9)
    aa=aa+1;
end
saveas(fig1,strcat('images/bursts CTAR/features comparison/frequencies_lessbins.png'));
saveas(fig1,strcat('images/bursts CTAR/features comparison/frequencies_lessbins.fig'));


%%
fig2=figure('units','normalized','outerposition',[0.01 0.01 0.99 0.99]);
sgtitle('Comparison fooof and bycycle frequencies')
subplot(3,2,1)
x=x_var_ecog;y=y_absdiff_ecog;
scatter(x,y,40,colors,'filled');%drawnow
alpha(0.5)
hold on
poly = polyfit(x,y,1);
y_fit = polyval(poly,x);
SST = sum((y - mean(y)).^2);
SSR = sum((y - y_fit).^2);
R2 = 1-(SSR/SST);
x_lim=xlim;
y_lim=ylim;
text( x_lim(1) + 0.3*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)),...
    sprintf('R2 = %.3f',R2),'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
plot(x,y_fit,'k-','LineWidth',2)
ylabel('abs diff ifreq - bycycle freq');title('ecog')
xlabel('Variance bycycle frequencies')

subplot(3,2,2)
x=x_var_dbs;y=y_absdiff_dbs;
scatter(x,y,40,[0.4,0.6,1],'filled');%drawnow
alpha(0.5)
hold on
poly = polyfit(x,y,1);
y_fit = polyval(poly,x);
SST = sum((y - mean(y)).^2);
SSR = sum((y - y_fit).^2);
R2 = 1-(SSR/SST);
x_lim=xlim;
y_lim=ylim;
text( x_lim(1) + 0.3*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)),...
    sprintf('R2 = %.3f',R2),'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
plot(x,y_fit,'k-','LineWidth',2)
ylabel('absolute difference ifreq - bycycle freq');title('dbs')
xlabel('Variance bycycle frequencies')

subplot(3,2,3)
x=x_var_ecog;y=y_zscore_ecog;
scatter(x,y,40,colors,'filled');%drawnow
alpha(0.5)
hold on
poly = polyfit(x,y,1);
y_fit = polyval(poly,x);
SST = sum((y - mean(y)).^2);
SSR = sum((y - y_fit).^2);
R2 = 1-(SSR/SST);
x_lim=xlim;
y_lim=ylim;
text( x_lim(1) + 0.3*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)),...
    sprintf('R2 = %.3f',R2),'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
plot(x,y_fit,'k-','LineWidth',2)
ylabel('zscore - bycycle freq')
xlabel('Variance bycycle frequencies')

subplot(3,2,4)
x=x_var_dbs;y=y_zscore_dbs;
scatter(x,y,40,[0.4,0.6,1],'filled');%drawnow
alpha(0.5)
hold on
poly = polyfit(x,y,1);
y_fit = polyval(poly,x);
SST = sum((y - mean(y)).^2);
SSR = sum((y - y_fit).^2);
R2 = 1-(SSR/SST);
x_lim=xlim;
y_lim=ylim;
text( x_lim(1) + 0.3*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)),...
    sprintf('R2 = %.3f',R2),'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
plot(x,y_fit,'k-','LineWidth',2)
ylabel('zscore - bycycle freq')
xlabel('Variance bycycle frequencies')

subplot(3,2,5)
x=x_f_ecog;y=y_b_ecog;
scatter(x,y,40,colors,'filled');%drawnow
alpha(0.5)
hold on
poly = polyfit(x,y,1);
y_fit = polyval(poly,x);
SST = sum((y - mean(y)).^2);
SSR = sum((y - y_fit).^2);
R2 = 1-(SSR/SST);
x_lim=xlim;
y_lim=ylim;
text( x_lim(1) + 0.3*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)),...
    sprintf('R2 = %.3f',R2),'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
plot(x,y_fit,'k-','LineWidth',2)
xlabel('Fooof frequencies')
ylabel('Bycycle freq')

subplot(3,2,6)
x=x_f_dbs;y=y_b_dbs;
scatter(x,y,40,[0.4,0.6,1],'filled');%drawnow
alpha(0.5)
hold on
poly = polyfit(x,y,1);
y_fit = polyval(poly,x);
SST = sum((y - mean(y)).^2);
SSR = sum((y - y_fit).^2);
R2 = 1-(SSR/SST);
x_lim=xlim;
y_lim=ylim;
text( x_lim(1) + 0.3*(x_lim(2)-x_lim(1)), y_lim(1) + 0.85*(y_lim(2)-y_lim(1)),...
    sprintf('R2 = %.3f',R2),'FontSize',10, 'FontWeight','bold','BackgroundColor', 'white', 'EdgeColor', 'black')
plot(x,y_fit,'k-','LineWidth',2)
xlabel('Fooof frequencies')
ylabel('Bycycle freq')

saveas(fig2,strcat('images/bursts CTAR/features comparison/trend frequencies.png'));
saveas(fig2,strcat('images/bursts CTAR/features comparison/trend frequencies.fig'));
close all

