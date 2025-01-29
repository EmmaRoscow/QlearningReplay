
function [h] = plot_reliability_error_over_sessions_x_replays(rat_names, ksamples, policy, cmap)

    colour_scheme;

    h = figure('Position', [680 558 680 420]); hold on

    for k = 1:length(ksamples)

        params.ksamples = ksamples(k);

        for iRat = 1:length(rat_names)

            rat = rat_names{iRat};

            % Load behavioural results
            load(['./data/behavioural_data/behaviour_' rat '.mat'])
            evalc(['behav_data = ' rat]);

            % Load optimised parameters
            load(['./data/q_learning_results/optimised_parameters/' rat '_' policy '.mat'])

            % Load reliability error per trial
            try
                load(['./data/q_learning_results/predictive_accuracy_' num2str(params.ksamples) '_replays.mat'])
                evalc(['reliability_error = predictive_accuracy.' policy '.' rat]);
            catch
                continue
            end

            % Normalise reliability error
            if params.ksamples == 0
                baseline(iRat) = mean(reliability_error);
            end
            reliability_error = reliability_error/baseline(iRat);

            % Group by training session
            for iTrial = 1:length(reliability_error)
                evalc(['session(iTrial) = find(iTrial > [0 cumsum(' rat '.n_trials-1)], 1, ''last'');']);
            end

            % Get normalised reliability error per session
            for iSession = 1:length(behav_data.n_trials)
                mean_reliability_error(k, iRat, iSession) = mean(reliability_error(session==iSession));
            end

            clear reliability_error session

        end

        % Plot mean +- standard error
        c = cmap(max([k*10, 1]), :);
        nSessions = length(behav_data.n_trials);
        mean_per_session = squeeze(mean(mean_reliability_error(k, :, :), 2));
        line([1:nSessions], mean_per_session,...
            'Color', c, 'LineWidth', 2)
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
    xlim([0.5 nSessions+0.5])
    set(gca, 'YTick', [0.1:0.1:0.5 1:5]);

    % Lines indicating learning stages
    yyaxis right; set(gca, 'YColor', 'none')
    line([15.5 15.5], [0 1], 'Color', colourscheme.learningstages)
    line([20.5 20.5], [0 1], 'Color', colourscheme.learningstages)
    drawnow

    % Text indicating learning stages
    opts = struct('HorizontalAlignment', 'center', 'Color', colourscheme.learningstages, 'FontSize', 14);
    t1 = text(8, 0.95, 'Initial learning', opts);
    t2 = text(18, 0.95, 'Revaluation', opts);
    t3 = text([21.5 21.5], [0.97 0.93], {'Rev-', 'ersal'}, opts);

    % Colorbar
    colormap(cmap);
    cbar = colorbar;
    cbar.Title.String = '# replays';
    cbar.Ticks = linspace(0, 1, length(ksamples));
    cbar.TickLabels = ksamples;
    cbar.Box = 'off';
    
end

