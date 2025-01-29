
colour_scheme;
axes_properties;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
params.seed = 0;
params.policy = 'none';
params.ksamples = 0;

plot_type = 'histogram';

for iRat = 1:length(rat_names)
    
    rat = rat_names{iRat};
    
    % Load behavioural results
    load(['./data/behavioural_data/behaviour_' rat '.mat'])
    evalc(['behav_data = ' rat]);
    
    % Load optimised parameters for no-replay baseline
    load(['./data/q_learning_results/optimised_parameters/' rat '_random.mat'])
    params.alpha = result{ksamples==0}(1);
    params.discount = result{ksamples==0}(2);
    params.epsilon = result{ksamples==0}(3);

    % Run DynaQ action prediction
    [predictive_accuracy, output] = dynaQ(behav_data, params, false);
    
    % Get session number for each trial
    evalc(['n_trials = [0 cumsum(' rat '.n_trials-1)]']);
    session_no{iRat} = NaN(1, n_trials(end));
    for iSession = 1:22
        session_no{iRat}(n_trials(iSession)+1 : n_trials(iSession+1)) = iSession;
    end
    session_no{iRat} = repmat(session_no{iRat}', 1, 3);
    
    % For each trial, get difference between predicted action probability
    % and overall observed action frequency
    action_prob_error{iRat} = output.action_prob - output.obs_bar(output.s_history, :);
    
end

all_action_prob_error = cell2mat(cellfun(@(x) x(:), action_prob_error, 'UniformOutput', false)');
all_session_no = cell2mat(cellfun(@(x) x(:), session_no, 'UniformOutput', false)');
    
% Bin predicted action errors into percentiles
[N, edges, bin] = histcounts(all_action_prob_error, -1:0.002:1);
    
% Plot distribution of residuals and the average session from which they
% come
figure('Position', [680 558 680 420]); hold on

switch plot_type

    case 'histogram'

        cmap = parula(2200);
        for iBin = 1:length(N)
            if N(iBin) > 0
                % Get average session
                mean_sess = mean(all_session_no(bin==iBin));
                % Plot
                bar(mean(edges(iBin:iBin+1)), N(iBin), 'FaceColor', cmap(round(mean_sess*100), :), 'EdgeColor', 'none', 'BarWidth', 0.004)
                blah(iBin) = mean_sess;
            end
        end

        % Add colorbar
        c = colorbar;
        c.Ticks = linspace(c.Limits(1), c.Limits(2), 22);
        c.Ticks = c.Ticks([1 5:5:20 22]);
        c.TickLabels = [1 5:5:20 22];
        c.Title.String = {'Average', 'session #'};

        % Formatting
        set(gca, 'YScale', 'log');
        xlim([-1 1])

        % Format y axis
        tick = get(gca, 'YTick');
        include = [1 , 10, 100, 1000];
        for i = 1:length(include)
            if ~ismember(include(i), tick); tick = [tick include(i)]; end
        end
        set(gca, 'YTick', sort(tick));
        set(gca, 'YTickLabel', sprintfc('%g', sort(tick)))

        % Axis labels
        xlabel({'Error between predicted action probability', 'and observed action frequency'});
        ylabel('Count (number of trials)')

        % Layout
        ax_pos = get(gca, 'Position');
        c.Position = [0.8373    0.2000    0.0392    0.6452];
        set(gca, 'Position', ax_pos);

    case 'scatter'
        
        for iSess = 1:22
            xjitter = linspace(iSess-0.25, iSess+0.25, sum(all_session_no==iSess));
            xjitter = xjitter(randperm(length(xjitter)));
            scatter(xjitter, all_action_prob_error(all_session_no==iSess), '.k')
        end

        % Formatting
        xlim([0 23])
        set(gca, 'XTick', [1 15 20 22])
        xlabel('Training session')
        ylabel({'Error between', 'predicted action probability', 'and observed action frequency'})
        
        % Lines indicating learning stages
        yyaxis right; set(gca, 'YColor', 'none')
        line([15.5 15.5], [0 1], 'Color', colourscheme.learningstages)
        line([20.5 20.5], [0 1], 'Color', colourscheme.learningstages)

        % Text indicating learning stages
        opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
        t1 = text(8, 0.95, 'Initial learning', opts);
        t2 = text(18, 0.95, 'Revaluation', opts);
        t3 = text([21.5 21.5], [0.97 0.93], {'Rev-', 'ersal'}, opts);
        
end

