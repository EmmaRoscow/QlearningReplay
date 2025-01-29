
colour_scheme;
axes_properties;
rat = 'Harry';
policies = {'random', 'expectedReward', 'RPE'};

% Load behavioural data
load(['./data/behavioural_data/behaviour_' rat '.mat'])
evalc(['behav_data = ' rat]);

% Load Q-values
load('./data/q_learning_results/Q_history_100_replays.mat')

for iPol = 1:length(policies)
    
    policy = policies{iPol};

    % Create figure
    h = figure('Position', [680 558 1360 420]); hold on

    % Plot baseline Q-values
    Q_history = Q_history_100_replays.none.Harry;
    h = plot_q_value_evolution(Q_history, behav_data, h);

    % Plot Q-values
    Q_history = Q_history_100_replays.(policy).Harry;
    h = plot_q_value_evolution(Q_history, behav_data, h);

    % Format axes and add lines for learning stages
    for i = 1:3
        subplot(1, 3, i)

        % If patches are obscuring lines, swap them
        ax = gca;
        if isa(ax.Children(2), 'matlab.graphics.primitive.Patch') && isa(ax.Children(end), 'matlab.graphics.chart.primitive.Line')
            ax.Children = flip(ax.Children);
        end
        clear ax

        % Axes
        n_trials = size(Q_history{1}, 1);
        xlim([0 n_trials+1])
        xlabel('Trial')
        set(gca, 'XTick', [1, n_trials])
        set(gca, 'YTick', [0 1])

        % Lines indicating learning stages
        line([1 1]*sum(behav_data.n_trials(1:15))-15+0.5, [0 1], 'Color', colourscheme.learningstages)
        line([1 1]*sum(behav_data.n_trials(1:20))-20+0.5, [0 1], 'Color', colourscheme.learningstages)

    end

    subplot(1, 3, 1); ylabel('Q-value')

    % Reduce height of plots to fit text above
    for i = 1:3
        subplot(1, 3, i); ax{i} = gca;
    end
    base = ax{1}.Position(2) + 0.05;
    height = ax{1}.Position(4) - 0.15;
    order = behav_data.reward_probs{1};
    order_idx = [find(strcmp(order, 'high')), find(strcmp(order, 'medium')), find(strcmp(order, 'low'))];
    for i = 1:3
        ax{i}.Position(4) = height;
        ax{i}.Position(2) = base;
        t = title(ax{i}, ['state = ' order{order_idx(i)}]);
        t.Position(2) = t.Position(2) + 0.1;
    end

    % Text indicating learning stages
    for i = 1:3
        opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
        initial_midpoint = (sum(behav_data.n_trials(1:15))-15+0.5)/2;
        revaluation_midpoint = (sum(behav_data.n_trials(1:20))-20+0.5)/2 + initial_midpoint;
        reversal_midpoint = (sum(behav_data.n_trials(1:22))-22+0.5)/2 + revaluation_midpoint - initial_midpoint;
        t1 = text(ax{i}, initial_midpoint, 1.05, 'Initial learning', opts);
        t2 = text(ax{i}, revaluation_midpoint, 1.05, 'Revaluation', opts);
        t3 = text(ax{i}, [reversal_midpoint reversal_midpoint], [1.07 1.03], {'Rev-', 'ersal'}, opts);
    end

    % Text indicating reward prob (state = high)
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.high, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{1}, initial_midpoint+100, 0.05, '0% reward', opts);
    text(ax{1}, revaluation_midpoint, 0.05, '0%', opts);
    text(ax{1}, reversal_midpoint, 0.05, '0%', opts);
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.medium, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{1}, initial_midpoint, 0.95, '50% reward', opts);
    text(ax{1}, revaluation_midpoint, 0.73, '50%', opts);
    text(ax{1}, reversal_midpoint, 0.75, '50%', opts);
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.low, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{1}, initial_midpoint, 0.5, '25% reward', opts);
    text(ax{1}, revaluation_midpoint, 0.63, '12.5%', opts);
    text(ax{1}, reversal_midpoint, 0.57, '87.5%', opts);

    % Text indicating reward prob (state = medium)
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.high, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{2}, initial_midpoint+100, 0.85, '75% reward', opts);
    text(ax{2}, revaluation_midpoint, 0.95, '87.5%', opts);
    text(ax{2}, reversal_midpoint, 0.8, '12.5%', opts);
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.medium, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{2}, initial_midpoint-100, 0.27, '0% reward', opts);
    text(ax{2}, revaluation_midpoint, 0.15, '0%', opts);
    text(ax{2}, reversal_midpoint, 0.15, '0%', opts);
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.low, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{2}, initial_midpoint-100, 0.5, '25% reward', opts);
    text(ax{2}, revaluation_midpoint, 0.4, '12.5%', opts);
    text(ax{2}, reversal_midpoint, 0.33, '87.5%', opts);

    % Text indicating reward prob (state = low)
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.high, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{3}, initial_midpoint, 0.67, '75% reward', opts);
    text(ax{3}, revaluation_midpoint, 0.5, '87.5%', opts);
    text(ax{3}, reversal_midpoint, 0.35, '12.5%', opts);
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.medium, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{3}, initial_midpoint, 0.6, '50% reward', opts);
    text(ax{3}, revaluation_midpoint, 0.4, '50%', opts);
    text(ax{3}, reversal_midpoint, 0.25, '50%', opts);
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.low, 'FontSize', 14, 'FontWeight', 'bold');
    text(ax{3}, initial_midpoint-100, 0.05, '0% reward', opts);
    text(ax{3}, revaluation_midpoint, 0.05, '0%', opts);
    text(ax{3}, reversal_midpoint, 0.05, '0%', opts);
    
end
