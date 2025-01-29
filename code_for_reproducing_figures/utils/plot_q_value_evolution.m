
function [h] = plot_q_value_evolution(Q_history, behav_data, h)

    colour_scheme;

    order = behav_data.reward_probs{1};
    order_idx = [find(strcmp(order, 'high')), find(strcmp(order, 'medium')), find(strcmp(order, 'low'))];
    % Iterate over states
    for i = 1:3
        figure(h); subplot(1, 3, i); hold on
        % Iterate over actions
        for j = 1:3
            % Lighter colours for little replay, darker colours for more
            % replay
            evalc(['c = colourscheme.' order{order_idx(j)}]);
            % If just an array, plot line
            if isa(Q_history, 'double')
                c = c/4*3;
                plot(squeeze(Q_history(:, order_idx(i), order_idx(j))), 'LineWidth', 3, 'Color', c)
            % If a cell array, plot shaded region
            else
                c = 1-((1-c)/2);
                % Get Q-value for current state-action pair
                Q = cellfun(@(x) x(:, order_idx(i), order_idx(j)), Q_history, 'UniformOutput', false);
                % Get max and min of distribution for every trial
                upper = max(cell2mat(Q), [], 2);
                lower = min(cell2mat(Q), [], 2);
                mean_q = mean(cell2mat(Q), 2);
                patch([1:length(upper') flip(1:length(upper'))], [upper' flip(lower')], 'k', 'EdgeColor', c, 'FaceColor', c, 'FaceAlpha', 0.5)
                % Add a line for the mean in case the distribution is so
                % thin that it's not visible
                plot(mean_q, 'LineWidth', 1.5, 'Color', c)
            end
            drawnow

        end
    end

end