
colour_scheme;
axes_properties;
rat = 'Harry';

% Load behavioural results
load(['./data/behavioural_data/behaviour_' rat '.mat'])
evalc(['behav_data = ' rat]);

% Create figure
h = figure('Position', [680 558 1360 420]); hold on

% Plot Q-values
order = behav_data.reward_probs{1};
order_idx = [find(strcmp(order, 'high')), find(strcmp(order, 'medium')), find(strcmp(order, 'low'))];
% Iterate over states
for i = 1:3
    figure(h); subplot(1, 3, i); hold on
    % Iterate over actions
    for j = 1:3

        evalc(['c = colourscheme.' order{order_idx(j)}]);

        % Do binomial test to get upper and lower boundaries of
        % proportion for the given action in the given state
        for iSession = 1:length(behav_data.actions)
            states = behav_data.actions{iSession}(1:end-1);
            actions = behav_data.actions{iSession}(2:end);
            n = sum(states == i);
            x = sum(states == i & actions == j);
            [~, pci] = binofit(x, n, 0.25);
            lower(iSession) = pci(1);
            upper(iSession) = pci(2);
        end

        % Plot 75% binomial confidence intervals
        patch([1:length(upper) flip(1:length(upper))], [upper flip(lower)], 'k', 'EdgeColor', c, 'FaceColor', c, 'FaceAlpha', 0.5)
        drawnow
            
    end
end


for i = 1:3
    subplot(1, 3, i)
        
    % Axes
    xlim([0.5 22.5])
    xlabel('Session')
    set(gca, 'XTick', [1, 22])
    set(gca, 'YTick', [0 1])
    
    % Lines indicating learning stages
    line([1 1]*15, [0 1], 'Color', colourscheme.learningstages)
    line([1 1]*20, [0 1], 'Color', colourscheme.learningstages)
    
end

subplot(1, 3, 1); ylabel('Frequency of action selection')

% Reduce height of plots to fit text above
for i = 1:3
    subplot(1, 3, i); ax{i} = gca;
end
base = ax{1}.Position(2) + 0.05;
height = ax{1}.Position(4) - 0.15;
for i = 1:3
    ax{i}.Position(4) = height;
    ax{i}.Position(2) = base;
    t = title(ax{i}, ['state = ' order{order_idx(i)}]);
    t.Position(2) = t.Position(2) + 0.1;
end

% Text indicating learning stages
for i = 1:3
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
    initial_midpoint = 7.5;
    revaluation_midpoint = 18;
    reversal_midpoint =21.5;
    t1 = text(ax{i}, initial_midpoint, 1.05, 'Initial learning', opts);
    t2 = text(ax{i}, [revaluation_midpoint revaluation_midpoint], [1.07 1.03], {'Reval-', 'uation'}, opts);
    t3 = text(ax{i}, [reversal_midpoint reversal_midpoint reversal_midpoint], [1.09 1.05 1.01], {'Re-', 'ver-', 'rsal'}, opts);
end

