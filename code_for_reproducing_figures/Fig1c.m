
colour_scheme;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};

for iRat = 1:length(rat_names)
    
    rat = rat_names{iRat};
    
    % Load behavioural results
    load(['./data/behavioural_data/behaviour_' rat '.mat'])
    
    % Get the arm identity chosen on each trial
    evalc(['actions = cellfun(@(x, y) arrayfun(@(z) x(z), y), ' rat '.reward_probs, ' rat '.actions, ''UniformOutput'', false)']);
    
    % Pool arm choices
    pooled_high{iRat} = cellfun(@(x) strcmp(x, 'high'), actions, 'UniformOutput', false);
    pooled_medium{iRat} = cellfun(@(x) strcmp(x, 'medium'), actions, 'UniformOutput', false);
    pooled_low{iRat} = cellfun(@(x) strcmp(x, 'low'), actions, 'UniformOutput', false);
    
    % Pool average arm choices
    mean_high(iRat, :) = cellfun(@(x) mean(strcmp(x, 'high')), actions);
    mean_medium(iRat, :) = cellfun(@(x) mean(strcmp(x, 'medium')), actions);
    mean_low(iRat, :) = cellfun(@(x) mean(strcmp(x, 'low')), actions);
    
end

figure('Position', [680 558 680 420]); hold on
    
% Line indicating chance at 33.3%
line([0 23], [1/3 1/3], 'Color', colourscheme.chance, 'LineStyle', '--', 'LineWidth', 3, 'HandleVisibility', 'off')

% Plot mean +- standard error
sem = std(mean_high) / sqrt(6);
errorbar(1:22, mean(mean_high), sem,...
    'Color', colourscheme.high, 'LineWidth', 2)
sem = std(mean_medium) / sqrt(6);
errorbar(1:22, mean(mean_medium), sem,...
    'Color', colourscheme.medium, 'LineWidth', 2)
sem = std(mean_low) / sqrt(6);
errorbar(1:22, mean(mean_low), sem,...
    'Color', colourscheme.low, 'LineWidth', 2)

% Perform chi-squared tests and plot significance
for iSession = 1:22
    n_trials = length(cell2mat(cellfun(@(x) x(iSession), pooled_high)));
    n_high = sum(cell2mat(cellfun(@(x) x(iSession), pooled_high)));
    n_medium = sum(cell2mat(cellfun(@(x) x(iSession), pooled_medium)));
    n_low = sum(cell2mat(cellfun(@(x) x(iSession), pooled_low)));
    [h, p(iSession), stats] = chi2gof([1 2 3], 'NBins', 3, 'Frequency', [n_high n_medium n_low], 'Expected', [n_trials/3 n_trials/3 n_trials/3]);
    if p(iSession) < 0.05
        text(iSession, 0.55, '*', 'FontSize', 24, 'HorizontalAlignment', 'center')
    end
end

% Formatting
axes_properties;
xlim([0.5 22.5])
ylim([0 0.7])

% Lines indicating learning stages
line([15.5 15.5], [0 1], 'Color', colourscheme.learningstages)
line([20.5 20.5], [0 1], 'Color', colourscheme.learningstages)

% Text indicating learning stages
opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
t1 = text(8, 0.65, 'Initial learning', opts);
t2 = text(18, 0.65, 'Revaluation', opts);
t3 = text([21.5 21.5], [0.67 0.63], {'Rev-', 'ersal'}, opts);

% Axis labels
xlabel('Training session'); set(gca, 'XTick', [1 15 20 22])
ylabel('Percentage of arm choices')
set(gca, 'YTick', [0 1/3 0.7]); set(gca, 'YTickLabel', {'0%', '33%', '70%'})

% Legend
l = legend('High prob. arm', 'Mid prob. arm', 'Low prob. arm');
l.Box = 'off'; l.Location = 'south'; l.Orientation = 'Horizontal';
