
colour_scheme;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
policies = {'noReplay', 'random', 'expectedReward', 'RPE', 'RPEProb'};
arm_order = {'high', 'medium', 'low'};
ksamples = 100;

% Load baseline predictive accuracy (no replay)
load('./data/q_learning_results/predictive_accuracy_0_replays.mat')
for iRat = 1:length(rat_names)
    rat = rat_names{iRat};
    evalc(['baseline(iRat) = mean(predictive_accuracy.noReplay.' rat ');']);
end

% Load reliability error per trial
load(['./data/q_learning_results/predictive_accuracy_' num2str(ksamples) '_replays.mat'])

figure('Position', [680 458 680 520]);

% Get change in reliability error for each policy
for iPol = 1:length(policies)
    
    policy = policies{iPol};
    pred_accuracy = cell(22, 3, 3);

    for iRat = 1:length(rat_names)
        
        rat = rat_names{iRat};
        rel = [];
        
        % Load behavioural results
        load(['./data/behavioural_data/behaviour_' rat '.mat'])
        evalc(['behav_data = ' rat]);

        % Get reliability error per trial for this policy
        evalc(['reliability_error = predictive_accuracy.' policy '.' rat]);

        % Normalise reliability error
        reliability_error = reliability_error/baseline(iRat);

        % Get the session number, state and action for each trial
        for iTrial = 1:length(reliability_error)
            evalc(['session(iTrial) = find(iTrial > [0 cumsum(' rat '.n_trials-1)], 1, ''last'');']);
        end
        actions = cell2mat(cellfun(@(x) x(2:end), behav_data.actions, 'UniformOutput', false));
        states = [0 actions(1:end-1)];
        
        % Get normalised reliability for each state-action pair
        for iSess = 1:length(behav_data.n_trials)
        
            for iState = 1:3
            
                for iAction = 1:3
                    
                    s = find(strcmp(behav_data.reward_probs{iSess}, arm_order{iState}));
                    a = find(strcmp(behav_data.reward_probs{iSess}, arm_order{iAction}));
        
                    % Get reliability error for all trials in corresponding
                    % session, state and action
                    rel = reliability_error(find(session==iSess & states==s & actions==a));
                        
                    % Add to something like pred_accuracy.policy.session.state.action to create
                    % a long vector
                    pred_accuracy{iSess, iState, iAction} = [pred_accuracy{iSess, iState, iAction} rel];
                    
                end
            end
        end
        clear session
    end
    
    % Set no-replay policy as baseline
    if strcmp(policy, 'noReplay')
        for iState = 1:3
            for iAction = 1:3
                no_replay(:, iState, iAction) = cellfun(@mean, pred_accuracy(:, iState, iAction));
            end
        end
        
    else
                    
        for iState = 1:3

            for iAction = 1:3

               % Plot the average over sessions
                subplot(3, 3, (iAction-1)*3 + iState); hold on
                evalc(['c = colourscheme.' policy ';']);
                to_plot =  cellfun(@mean, pred_accuracy(:, iState, iAction));
                to_plot = to_plot ./ no_replay(:, iState, iAction);
                xdata = find(~isnan(to_plot));
                ydata = to_plot(~isnan(to_plot));
                plot(xdata, ydata, 'color', c, 'LineWidth', 2.5)
                drawnow

            end
        end
    end
end

% Formatting etc.
axes_properties;
xlim([0.5 22.5])
for iState = 1:3
    for iAction = 1:3
        subplot(3, 3, (iState-1)*3 + iAction); hold on
        axis tight
        tick = get(gca, 'YTick'); set(gca, 'YTick', [0 tick(end)]);
        tick = get(gca, 'YTickLabel');
        set(gca, 'YTickLabel', cellfun(@(x) [num2str(str2num(x)*100) '%'], tick, 'UniformOutput', false));
        set(gca, 'XTick', [1 15 22])
        if iState < 3; set(gca, 'XTickLabel', []); end
        x = xlim; line(x, [1 1], 'LineStyle', '--', 'Color', colourscheme.chance)
    end
end

% Titles to indicate action
subplot(3, 3, 1); title({'Action = ', 'high prob. arm'}, 'fontsize', 16)
subplot(3, 3, 2); title({'Action = ', 'mid prob. arm'}, 'fontsize', 16)
subplot(3, 3, 3); title({'Action = ', 'low prob. arm'}, 'fontsize', 16)

% Y labels to indicate action
subplot(3, 3, 1); ylabel({'State = ', 'high prob.', 'arm'}, 'fontsize', 16, 'fontweight', 'bold')
subplot(3, 3, 4); ylabel({'State = ', 'mid prob.', 'arm'}, 'fontsize', 16, 'fontweight', 'bold')
subplot(3, 3, 7); ylabel({'State = ', 'low prob.', 'arm'}, 'fontsize', 16, 'fontweight', 'bold')

subplot(3, 3, 8); xlabel('Training session')

% Reduce spaces between subplots
drawnow
subplot(3, 3, 1); ax{1} = gca; height = ax{1}.Position(4) + 0.05; width = ax{1}.Position(3) + 0.03;
for i = 2:9
    subplot(3, 3, i)
    ax{i} = gca;
end
for i = 1:3
    ax{i}.Position(2) = ax{i}.Position(2) - 0.05   ;
    ax{i}.Position(3) = width;
    ax{i}.Position(4) = height;
end
drawnow
for i = 4:6
    ax{i}.Position(2) = ax{i}.Position(2) - 0.025;
    ax{i}.Position(3) = width;
    ax{i}.Position(4) = height;
end
drawnow
for i = 7:9
    ax{i}.Position(3) = width;
    ax{i}.Position(4) = height;
end
drawnow
for i = [3 6 9]
    ax{i}.Position(1) = ax{i}.Position(1) + 0.02;
end
for i = [2 5 8]
    ax{i}.Position(1) = ax{i}.Position(1) + 0.04;
end
for i = [1 4 7]
    ax{i}.Position(1) = ax{i}.Position(1) + 0.06;
    ax{i}.YLabel.Position(1) = ax{i}.YLabel.Position(1)-7;
end

% Y label
text(ax{4}, -15, 1, 'Error score relative to no-replay baseline', 'rotation', 90, 'horizontalalignment', 'center',...
    'FontName', ax{8}.XLabel.FontName,...
    'FontSize', ax{8}.XLabel.FontSize,...
    'Color', ax{8}.XLabel.Color);
