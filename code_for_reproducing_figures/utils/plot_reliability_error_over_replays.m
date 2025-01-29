
function [h, pval, ksamples] = plot_reliability_error_over_replays(rat_names, policies, shuffled)

    colour_scheme;
    pred_accuracy_combined = cell(1, length(policies));

    h = figure('Position', [680 558 680 420]); hold on

    for iPol = 1:length(policies)

        policy = policies{iPol};

        pred_error = [];
        pred_accuracy_combined{iPol} = [];

        for iRat = 1:length(rat_names)

            rat = rat_names{iRat};

            % Load optimised parameters
            if shuffled
                load(['./data/q_learning_results/optimised_parameters/shuffled/' rat '_' policy '.mat'])
            else
                load(['./data/q_learning_results/optimised_parameters/' rat '_' policy '.mat'])
            end

            % Normalise
            baseline = test_score(ksamples==0);
            test_score = test_score/baseline;

            % Store
            pred_accuracy_combined{iPol}(:, iRat) = test_score';

            % Load reliability error for every trial
%             for k = 1:length(ksamples)
% 
%                 % Load reliability error for every trial
%                 if shuffled
%                     load(['./data/shuffled/predictive_accuracy_' num2str(ksamples(k)) '_replays.mat'])
%                 else
%                     load(['./data/predictive_accuracy_' num2str(ksamples(k)) '_replays.mat'])
%                 end
% 
%                 % Get normalised predictive accuracy
%                 evalc(['pred_accuracy_combined{iPol}{k}{iRat} = predictive_accuracy.' policy '.' rat '/baseline']);
% 
%             end

        end
        
        
        % Calculate significance using linear mixed-effects model
        for k = 1:length(ksamples) 
        
            if ksamples(k) == 0

                pval(iPol, k) = 1;
                intercept(iPol, k) = 1;
                sem(iPol, k) = 0;

            else

                % Concatenate predictive accuracy for given replay policy and number of replays, for all rats, and subtract 1 to test if sigifnicantly different from 0
                pred_acc_replay = pred_accuracy_combined{iPol}(k, :)' - 1;
                % Create vector of rat IDs as additional variable
                rat_id = cell2mat(cellfun(@(x, y) ones(length(x), 1) * y, mat2cell(pred_accuracy_combined{1}(1, :)', ones(length(pred_accuracy_combined{1}(1, :)), 1), 1), mat2cell([1:length(rat_names)]', ones(length(rat_names), 1), 1), 'UniformOutput', false));
                % Add variables to table
                tbl = table(pred_acc_replay, rat_id);
                % Fit linear mixed-effects model
                lme = fitlme(tbl, 'pred_acc_replay ~ 1 + (1|rat_id)');
                % Get one-sided p-value for fixed effect of replay
                idx = find(cellfun(@(x) strcmp(x, '(Intercept)'), lme.Coefficients.Name));
                if lme.Coefficients.Estimate(idx) < 0
                    pval(iPol, k) = lme.Coefficients.pValue( idx ) / 2;
                else
                    pval(iPol, k) = 1;
                end
                % Get the estimated mean
                intercept(iPol, k) = lme.Coefficients.Estimate(idx) + 1;
                % Get standard error of the mean
                sem(iPol, k) = lme.Coefficients.SE(idx);
                
            end
            
        end

        % Plot mean +- standard error
        evalc(['c = colourscheme.' policy]);
        errorbar(ksamples + -0.75+0.5*iPol, intercept(iPol, :), sem(iPol, :),...
            'Color', c, 'LineWidth', 2)
        
    end
        
    % Formatting
    axes_properties;

    % Axes
    xlim([-3 ksamples(end)+3])
    y = ylim; y(1) = min([y(1) 0.5]); y(2) = max([y(2) 1.5]); set(gca, 'YLim', y)
    tick = get(gca, 'YTick'); set(gca, 'YTick', [tick(1) 1 tick(end)])
    ylabel('Normalised error score')
    xlabel('Number samples replayed')
    set(gca, 'Position', [0.26 0.144 0.645 0.781])
    y = ylim;
    ax = gca; ax.YAxis.Label.Position = [-30 mean(y) -1];
    
    % Y-axis annotations
    txt{1} = text(-15, mean([mean(y) y(2)]), {'Less accurate', 'than no replay'});
    txt{2} = text(-15, mean([mean(y) y(1)]), {'More accurate', 'than no replay'});
    for t = 1:2
        txt{t}.FontName = 'Arial';
        txt{t}.FontSize = 12;
        txt{t}.Rotation = 90;
        txt{t}.HorizontalAlignment = 'center';
    end

    % Line to indicate 100%
    x = xlim; line(x, [1 1], 'LineStyle', '--', 'Color', colourscheme.chance)
        
    % Legend
    label_mapping = {{'noReplay', 'random', 'expectedReward', 'RPE', 'RPEProb'},...
        {'No replay', 'Random', 'Reward-biased', 'RPE-prioritised', 'RPE-proportional'}};
    for iPol = 1:length(policies)
        legend_label(iPol) = label_mapping{2}(strcmp(label_mapping{1}, policies{iPol}));
    end
    l = legend(legend_label);
    l.Location = 'northwest';
    l.FontSize = 12;
    l.Title.String = 'Replay policy';
    l.Box = 'off';
    
end
