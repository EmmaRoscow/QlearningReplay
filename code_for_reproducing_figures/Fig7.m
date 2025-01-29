
rat = 'Harry';

% Load behavioural results
load(['./data/behavioural_data/behaviour_' rat '.mat'])
evalc(['behav_data = ' rat]);

% Load optimised parameters for no-replay baseline
load(['./data/q_learning_results/optimised_parameters/' rat '_RPE.mat'])

% Set parameters
params.ksamples = 20;
params.policy = 'RPE';
optimised_params = result{ksamples==params.ksamples};
params.alpha = optimised_params(1);
params.discount = optimised_params(2);
params.epsilon = optimised_params(3);
if length(optimised_params) >=4
    params.recency = optimised_params(4);
end
if length(optimised_params)==5
    params.rpe_recency = optimised_params(5);
end

% Run dynaQ
[~, output] = dynaQ(behav_data, params, false);

replay_history_abs = output.replay_history_abs;
n_trials = behav_data.n_trials;

% Plot
fig = figure('Position', [680 558 680 420]);
    
% Plot the trial from which sample is taken
subplot(1, 5, 1:4); hold on
for i = 1:length(replay_history_abs)
    xjitter = linspace(-0.1, 0.1, length(replay_history_abs{i}))+i; xjitter = xjitter(randperm(length(xjitter)));
    scatter(xjitter, replay_history_abs{i}, 3, 'k')
end

% Plot limit (most recent trial)
stairs(0.5:1:length(n_trials)+0.5, [cumsum(n_trials-1) sum(n_trials-1)], 'Color', [0.5 0.5 0.5])

% Axes
xlabel({'Replay following session number', ' '})
ylabel({'Trial from which', 'replayed sample is taken'})
xlim([0.5 length(n_trials)+0.5])
ylim([0 sum(n_trials-1)+1])
drawnow

% Plot number of times each trial is replayed
subplot(1, 5, 5);
histogram(cell2mat(replay_history_abs'), 1:sum(n_trials-1), 'Orientation', 'horizontal', 'FaceColor', 'k', 'EdgeColor', 'none')
set(gca, 'YTickLabel', [])
xlabel({'Total number', 'times replayed'})
axis tight; ylim([0 sum(n_trials-1)+1])
box off
drawnow

% Align plots
fig.Children(2).Position(1) = fig.Children(2).Position(1) + 0.05;
fig.Children(2).Position(3) = fig.Children(2).Position(3) - 0.05;
fig.Children(2).Position(2) = fig.Children(1).Position(2);
fig.Children(2).Position(4) = fig.Children(1).Position(4);

