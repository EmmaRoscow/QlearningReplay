
colour_scheme;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
mean_optimality = [];
pooled_ntrials = [];
pooled_noptimal = [];

% Load behavioural results from all rats
for iRat = 1:length(rat_names)
    
    % Load behavioural results
    rat = rat_names{iRat};
    load(['./data/behavioural_data/behaviour_' rat '.mat'])
    evalc(['behav_data = ' rat]);
    
    % Pool results
    mean_optimality(iRat, :) =  behav_data.proportion_optimal;
    if isempty(pooled_ntrials)
        pooled_ntrials = behav_data.n_trials;
    else
        pooled_ntrials = pooled_ntrials + behav_data.n_trials;
    end
    if isempty(pooled_noptimal)
        pooled_noptimal = behav_data.n_trials .* behav_data.proportion_optimal;
    else
        pooled_noptimal = pooled_noptimal + behav_data.n_trials .* behav_data.proportion_optimal;
    end
    
    
end

% One-proportion z-test of difference from 50% chance level
p_0 = 0.5;
p = pooled_noptimal ./ pooled_ntrials;
n = pooled_ntrials;
z = (p-p_0) ./ (sqrt(p_0.*(1-p_0)./n));
p = 1-normcdf(z);

% Plot mean +- standard error
figure('Position', [680 558 680 420]); hold on
sem = std(mean_optimality) / sqrt(length(rat_names));
errorbar(1:22, mean(mean_optimality), sem,...
    'k', 'LineWidth', 2)

% Plot z-test significance with Bonferroni correction
for iSession = 1:22
    if p(iSession) < 0.05/length(pooled_ntrials)
        text(iSession, 0.8, '*', 'FontSize', 24, 'HorizontalAlignment', 'center')
    end
end

% Formatting
axes_properties;
xlim([0.5 22.5])
ylim([0 1])

% Lines indicating chance at 33.3% and 50%
opts = struct('Color', colourscheme.chance, 'LineStyle', '--', 'LineWidth', 3);
line([0 23], [1/3 1/3], opts)
line([0 23], [0.5 0.5], opts)

% Lines indicating learning stages
line([15.5 15.5], [0 1], 'Color', colourscheme.learningstages)
line([20.5 20.5], [0 1], 'Color', colourscheme.learningstages)

% Text indicating learning stages
opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
t1 = text(8, 0.95, 'Initial learning', opts);
t2 = text(18, 0.95, 'Revaluation', opts);
t3 = text([21.5 21.5], [0.97 0.93], {'Rev-', 'ersal'}, opts);

% Axis labels
xlabel('Training session'); set(gca, 'XTick', [1 15 20 22])
ylabel('Percentage optimal choices')
set(gca, 'YTick', [0 1/3 1/2 1]); set(gca, 'YTickLabel', {'0%', '33%', '50%', '100%'})
