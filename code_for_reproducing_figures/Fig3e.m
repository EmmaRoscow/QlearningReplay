
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};

% Load results of perturbation analysis
load('./data/q_learning_results/perturbation_analysis_results.mat')

% Plot
figure('Position', [680 558 450 420]); hold on
imagesc((perturbation_factor-1)*100, 1:3*length(rat_names), [perturbed_reliability.alpha; perturbed_reliability.discount; perturbed_reliability.epsilon]);
set(gca, 'YDir', 'reverse')
axis tight

% Formatting
axes_properties;

% Lines between parameters
x = xlim; line(x, [1 1] * length(rat_names) + 0.5, 'LineWidth', 4, 'Color', 'w')
x = xlim; line(x, [2 2] * length(rat_names) + 0.5, 'LineWidth', 4, 'Color', 'w')

% Axes
set(gca, 'YTick', mean([1, length(rat_names)]) + length(rat_names)*[0 1 2])
set(gca, 'YTickLabel', {'\alpha', '\gamma', '\epsilon'})
ylabel('Parameter')
ax = gca; ax.YAxis.FontSize = 36; ax.YAxis.Label.FontSize = 16;
ax.XAxis.TickLabelFormat = '%g%%';
xlabel('Perturbation from optimised value')
 y = ylim; yyaxis right; ylim(y); set(gca, 'YDir', 'reverse', 'YColor', 'k')
 set(gca, 'YTick', 1:3*length(rat_names))
set(gca, 'YTickLabel', repmat(cellfun(@(x) x(1), rat_names), 1, 3)')
ylabel('Rat')

% Colorbar
c = colorbar;
c.Label.String = {'Normalised', 'error score'};
