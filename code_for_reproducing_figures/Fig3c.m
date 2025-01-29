
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
params.seed = 0;
params.policy = 'none';
params.ksamples = 0;

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
    [reliability_error{iRat}] = dynaQ(behav_data, params, false);
    
end

% Plot distribution of reliability errors
figure('Position', [680 558 230 420])
boxplot(cell2mat(reliability_error)', cell2mat(cellfun(@(x, y) ones(length(x), 1)*y, reliability_error, num2cell(1:6), 'UniformOutput', false)'),...
    'PlotStyle', 'compact', 'colors', 'k', 'Whisker', inf);

% Formatting
axes_properties;
box off;

% Format axes
set(gca, 'YScale', 'log')
tick = get(gca, 'YTick');
include = [0.1, 1 , 10];
for i = 1:length(include)
    if ~ismember(include(i), tick); tick = [tick include(i)]; end
end
set(gca, 'YTick', sort(tick));
set(gca, 'YTickLabel', sprintfc('%g', sort(tick)))
set(gca, 'XTick', 1:length(rat_names)); set(gca, 'XTickLabel', cellfun(@(x) x(1), rat_names)');
xlabel('Rat'); ax = gca; ax.XLabel.Position(2) = ax.XLabel.Position(2) - 10;
ylabel('Error score');
