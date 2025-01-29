
function [h] = plot_reliability_error_over_sessions(rat_names, policies, shuffled, k)

    colour_scheme;
    
    h = figure('Position', [680 558 680 420]); hold on

    for iPol = 1:length(policies)

        policy = policies{iPol};

        for iRat = 1:length(rat_names)

            rat = rat_names{iRat};

            % Load behavioural results
            load(['./data/behavioural_data/behaviour_' rat '.mat'])
            evalc(['behav_data = ' rat]);

            % Load reliability error per trial
            if shuffled
                load(['./data/q_learning_results/shuffled/predictive_accuracy_' num2str(k) '_replays.mat'])
            else
                load(['./data/q_learning_results/predictive_accuracy_' num2str(k) '_replays.mat'])
            end
            evalc(['reliability_error = predictive_accuracy.' policy '.' rat]);

            % Normalise reliability error
            if strcmp(policy, 'noReplay')
                baseline(iRat) = mean(reliability_error);
            end
            reliability_error = reliability_error/baseline(iRat);

            % Group by training session
            for iTrial = 1:length(reliability_error)
                evalc(['session(iTrial) = find(iTrial > [0 cumsum(' rat '.n_trials-1)], 1, ''last'');']);
            end

            % Get normalised reliability error per session
            for iSession = 1:length(behav_data.n_trials)
                mean_reliability_error(iPol, iRat, iSession) = mean(reliability_error(session==iSession));
            end

            clear reliability_error session

        end

        % Plot mean +- standard error
        evalc(['c = colourscheme.' policy]);
        nSessions = length(behav_data.n_trials);
        sem_per_session = squeeze(std(mean_reliability_error(iPol, :, :), [], 2)) / sqrt(length(rat_names));
        mean_per_session = squeeze(mean(mean_reliability_error(iPol, :, :), 2));
        line([1:nSessions], mean_per_session,...
            'Color', c, 'LineWidth', 2)
        errorbar([1:nSessions] + -0.3+0.1*iPol, mean_per_session, sem_per_session,...
           '.',  'Color', c, 'LineWidth', 0.75, 'HandleVisibility', 'off')
        drawnow

    end

    % Formatting
    set(0, 'DefaultAxesFontSize', 16)
    set(0, 'DefaultAxesLineWidth', 1.5)
    set(0, 'DefaultAxesFontName', 'Arial')

    % Format axes
    set(gca, 'YScale', 'log')
    xlabel('Training session')
    ylabel('Normalised error score');
    xlim([0 nSessions+1])
    set(gca, 'YTick', [0.1:0.1:0.5 1:5]);

    % Lines indicating learning stages
    yyaxis right; set(gca, 'YColor', 'none')
    line([15.5 15.5], [0 1], 'Color', colourscheme.learningstages)
    line([20.5 20.5], [0 1], 'Color', colourscheme.learningstages)

    % Text indicating learning stages
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
    t1 = text(8, 0.95, 'Initial learning', opts);
    t2 = text(18, 0.95, 'Revaluation', opts);
    t3 = text([21.5 21.5], [0.97 0.93], {'Rev-', 'ersal'}, opts);

    % Title
    if k == 1 s = 'sample';
    else s = 'samples'; end
    title(sprintf('%i %s replayed between sessions', k, s), 'FontWeight', 'bold')

    % Legend
    label_mapping = {{'noReplay', 'random', 'expectedReward', 'RPE', 'RPEProb'},...
        {'No replay', 'Random', 'Reward-biased', 'RPE-prioritised', 'RPE-proportional'}};
    for iPol = 1:length(policies)
        legend_label(iPol) = label_mapping{2}(strcmp(label_mapping{1}, policies{iPol}));
    end
    l = legend(legend_label);
    l.Location = 'southwest';
    l.Title.String = 'Replay policy';
    l.Box = 'off';
    
end
