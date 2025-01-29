
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
params.seed = 0;
params.policy = 'none';
params.ksamples = 0;
colour_scheme;

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
    
    % Noramlise reliability error
    reliability_error{iRat} = reliability_error{iRat}/mean(reliability_error{iRat});
    
    % Group by training session
    for iTrial = 1:length(reliability_error{iRat})
        evalc(['session{iRat}(iTrial) = find(iTrial > [0 cumsum(' rat '.n_trials-1)], 1, ''last'');']);
    end
    
end

% Plot reliability errors
figure('Position', [680 558 680 420]); hold on
for iRat = 1:length(rat_names)
    xdata = session{iRat};
    xjitter = rand(1, length(xdata))*0.5 - 0.25;
    xdata = xdata + xjitter;
    scatter(xdata, reliability_error{iRat},...
        5, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceAlpha', 0.25)
end

% Formatting
axes_properties;
    
% Plot mean and standard error
for iSession = 1:22
    session_rel = cell2mat(cellfun(@(x, y) x(y==iSession), reliability_error, session, 'UniformOutput', false));
    session_mean(iSession) = mean(session_rel);
    session_sem(iSession) = std(session_rel) / sqrt(length(session_rel));
end
errorbar(session_mean, session_sem,...
    'k', 'LineWidth', 2)

% Format axes
set(gca, 'YScale', 'log')
tick = get(gca, 'YTick');
include = [0.1, 1 , 10];
for i = 1:length(include)
    if ~ismember(include(i), tick); tick = [tick include(i)]; end
end
set(gca, 'YTick', sort(tick));
set(gca, 'YTickLabel', sprintfc('%g', sort(tick)))
xlabel('Training session')
ylabel('Normalised error score');
xlim([0 23])
set(gca, 'XTick', [1 15 20 22])

% Add space to top to fit text
y = ylim;
y(2) = y(2)*5;
ylim(y);

% Lines indicating learning stages
yyaxis right; set(gca, 'YColor', 'none')
line([15.5 15.5], [0 1], 'Color', colourscheme.learningstages)
line([20.5 20.5], [0 1], 'Color', colourscheme.learningstages)

% Text indicating learning stages
opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
t1 = text(8, 0.95, 'Initial learning', opts);
t2 = text(18, 0.95, 'Revaluation', opts);
t3 = text([21.5 21.5], [0.97 0.93], {'Rev-', 'ersal'}, opts);
