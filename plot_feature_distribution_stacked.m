close all; clear all; clc;
%%
PATH_DATA = 'Z:\DBS';

DATE = datestr(now, 'yyyymmdd');
format long

% Scegli i dati da analizzare
n_sub_PD_DBS = [3003,3006,3008,3010:3012,3014,3015,3018,3020:3022,3024,3025,3027,3028];
n_sub_PD_DBS = arrayfun(@(x) sprintf('%04d', x), n_sub_PD_DBS, 'UniformOutput', false);
SUBJECTS = n_sub_PD_DBS;

tab_areas = readtable("HCPMMP1toAreas.txt");

%% Analisi
measures = {'duration','frequency','n_bursts','perc_bursts','volt_amp','time_rdsym','time_ptsym',};
measures_toplot = {'duration','frequency','n_bursts','probability','volt_amp','time_rdsym','time_ptsym',};
areas = unique(tab_areas.area);
areas(11) = {'dbs'};

% Ciclo sulle aree
for a = 1:numel(areas)
    area = areas{a};
    disp(strcat('-- ', area))
    if strcmp(area, 'LTC') || strcmp(area, 'IFOC'); continue; end

    % Inizializza matrici per salvare percentuali
    perc_positive_significant = zeros(1, numel(measures));
    perc_negative_significant = zeros(1, numel(measures));
    perc_not_significant = zeros(1, numel(measures));
    
    for m = 1:numel(measures)
        measure = measures{m};
        disp(strcat('-- Measure: ', measure))

        % -----------------PUT HERE WHAT TO PLOT!!-------------------------
        comparison = 'rebound'; % compare vs baseline
        stat_toEvaluate = 'ttest_t_rebound';
        pv_toEvaluate = 'ttest_pv_rebound';
        % -----------------------------------------------------------------

        %% Get data 
        total_electrodes = 0;
        num_positive_significant = 0;
        num_negative_significant = 0;
        num_not_significant = 0;

        ii = 1:numel(SUBJECTS);
        for i = ii
            if i == 6; continue; end
            SUBJECT = strcat('DBS', string(SUBJECTS(i)));
            disp(strcat('Now running i= ', string(i), '   aka: ', SUBJECT))

            PATH_ANNOT = strcat(PATH_DATA, filesep, SUBJECT, filesep, 'Preprocessed data\Sync\annot');
            electrode = bml_annot_read(strcat(PATH_ANNOT, filesep, SUBJECT, '_electrode'));
            electrode = electrode(:, {'id', 'starts', 'ends', 'duration', 'electrode', 'connector', 'port', 'HCPMMP1_label_1', 'HCPMMP1_weight_1'});
            cfg = [];
            cfg.decodingtype = 'basic';
            electrode = bml_getEcogArea(cfg, electrode);

            % Get BURSTS DATA
            tab_stats = readtable(strcat('annot/general CTAR/areas/', SUBJECT, " ", 'bycycle features comparison.txt'));

            if strcmp(area, 'dbs')
                tab_stats_area = tab_stats( (startsWith(tab_stats.label,'dbs') & strcmp(tab_stats.measure, measure)) ,:);
            else
                area_channels = electrode.electrode(strcmp(electrode.HCPMMP1_area, area)); 
                tab_stats_area = tab_stats( (ismember(tab_stats.label, area_channels) & strcmp(tab_stats.measure, measure)) ,:);
            end

            % Conta il numero di elettrodi analizzati per l'area
            electrodes_analysis = unique(tab_stats_area.label);
            total_electrodes = total_electrodes + numel(electrodes_analysis);

            % Conta i valori significativi positivi, negativi e non significativi
            for e = 1:numel(electrodes_analysis)
                tab_stats_electrode = tab_stats_area(strcmp(tab_stats_area.label, electrodes_analysis{e}),:);

                stat_comparison = tab_stats_electrode.(stat_toEvaluate);
                pv_comparison = tab_stats_electrode.(pv_toEvaluate);

                if pv_comparison < 0.05
                    if stat_comparison > 0
                        num_positive_significant = num_positive_significant + 1;
                    else
                        num_negative_significant = num_negative_significant + 1;
                    end
                else
                    num_not_significant = num_not_significant + 1;
                end
            end
        end

        % Calcola le percentuali
        perc_positive_significant(m) = (num_positive_significant / total_electrodes) * 100;
        perc_negative_significant(m) = (num_negative_significant / total_electrodes) * 100;
        perc_not_significant(m) = (num_not_significant / total_electrodes) * 100;
    end

    % Crea la figura per l'area
    fig = figure('units', 'normalized', 'outerposition', [0.03 0.03 0.7 0.5]);
    hold on
    %sgtitle(strcat('Proportion of significant results in ', area), 'FontWeight', 'bold')

    % Matrice per barre impilate
    stacked_data = [perc_not_significant; perc_positive_significant; perc_negative_significant]';

    % Crea le barre impilate
    h = bar(stacked_data, 'stacked');

    set(h(1), 'FaceColor', [0.9 0.9 0.9]);  % Non significativi - grigio chiaro
    set(h(2), 'FaceColor', [0.8 0.1 0.2]);        % Positivi significativi - rosso
    set(h(3), 'FaceColor', [0.2 0.2 0.8]);

    % Imposta colori (bianco per non significativi, rosso per positivi, blu per negativi)
    %colormap([0.9 0.9 0.9; 1 0 0; 0 0 1])

    % Assegnazione delle etichette
    set(gca, 'xtick', 1:numel(measures));
    set(gca, 'xticklabel', strrep(measures_toplot,'_',' '),'FontSize',16,'XTickLabelRotation', 30)%, 'XTickLabelRotation', 30)
    ylabel('Percentage (%)');
    ylim([0 100]);
    legend({'Not significant', 'Positive t value', 'Negative t value'}, 'Location', 'northeastoutside')

    % Salva la figura
    saveas(fig, strcat('images/bursts CTAR/areas/group/Proportion_significant_results',area,' ',comparison,'.png'))
    saveas(fig, strcat('images/bursts CTAR/areas/group/Proportion_significant_results',area,' ',comparison,'.fig'))
    close(fig)
end
