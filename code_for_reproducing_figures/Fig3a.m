
colour_scheme;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
params.seed = 0;
params.policy = 'none';
params.ksamples = 0;

% Initialise variables
mean_predicted_prob = cell(1, length(rat_names));
mean_observed_prob = cell(1, length(rat_names));
mean_session_no = cell(1, length(rat_names));


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
    session_no = NaN(1, n_trials(end));
    for iSession = 1:22
        session_no(n_trials(iSession)+1 : n_trials(iSession+1)) = iSession;
    end
    
    % Bin predicted action probabilities into percentiles
    edges = prctile(output.action_prob(:), [0:100]);
    [~, ~, bin] = histcounts(output.action_prob, edges);
    
    for iBin = 1:100
        
        % Get average of predicted probabilities
        mean_predicted_prob{iRat}(iBin) = mean(output.action_prob(bin==iBin));
        
        % Get average of actual probabilities
        [row, col] = find(bin==iBin);
        mean_observed_prob{iRat}(iBin) = mean(output.a_history(row) == col');
        
        % Get average session
        mean_session_no{iRat}(iBin) = mean(session_no(row));
        
    end
    
    % Get per-rat correlation of reliability diagram
    [rho(iRat), p(iRat)] = corr(mean_predicted_prob{iRat}', mean_observed_prob{iRat}');
    fprintf('Pearson''s correlation of %.2f between predicted and observed probabilities for rat %s\n',...
    rho(iRat), rat);

end

% Create vector of rat IDs
rat_id = cell2mat(cellfun(@(x, y) ones(length(x), 1) * y, mean_predicted_prob', mat2cell([1:length(rat_names)]', ones(length(rat_names), 1), 1), 'UniformOutput', false));
           
% Pool data
mean_predicted_prob = cell2mat(mean_predicted_prob)';
mean_observed_prob = cell2mat(mean_observed_prob)';
mean_session_no = cell2mat(mean_session_no)';

% Perform mixed-effects model
tbl = table(mean_observed_prob, mean_predicted_prob, rat_id);

% Fit linear mixed-effects model
lme = fitlme(tbl, 'mean_observed_prob ~ mean_predicted_prob-1 + (mean_predicted_prob-1|rat_id)');

% Get R-squared
rsquared = lme.Rsquared.Ordinary;

% Plot binned predicted and observed action probabilities
figure('Position', [680 558 680 420]); hold on
scatter(mean_predicted_prob, mean_observed_prob, 'k', 'HandleVisibility', 'off')

% Plot regression line
plot([0 1], [0 1]*lme.Coefficients.Estimate,...
    'k', 'LineWidth', 2)

% Plot perfect correlation
plot([0 1], [0 1], 'Color', colourscheme.chance, 'LineWidth', 2, 'LineStyle', '--')

% Legend
leg_text = sprintf('Regression (R^2 = %.2f)', rsquared');
l = legend(leg_text, 'Perfect correlation');
l.Box = 'off'; l.Location = 'southeast';
l.FontSize = 14;

% Formatting
axes_properties;
xlim([0 1])
ylim([0 1])

% Axis labels
xlabel('Predicted action probability'); set(gca, 'XTick', [0 1])
ylabel('Observed action probability'); set(gca, 'YTick', [0 1]); 
