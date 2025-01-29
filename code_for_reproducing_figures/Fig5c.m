
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
policies = {'noReplay', 'random', 'expectedReward', 'RPE', 'RPEProb'};
shuffled = false;
n_replays = 100;

% Plot
[h] = plot_reliability_error_over_sessions(rat_names, policies, shuffled, n_replays);

% Adjust size
h.Position(3) = h.Position(3)/3*2;

% Adjust text
hText = findobj(gcf, 'Type', 'Text');
hText(3).String = 'Reval-';
hText(3).Position(2) = 0.97;
hText(5) = copyobj(hText(3), gca);
hText(5).String = 'uation';
hText(5).Position(2) = 0.93;
hText(2).String = 'Re-';
hText(1).Position(1) = hText(1).Position(1)+0.2;
hText(2).Position(1) = hText(2).Position(1)+0.2;
hText(2).Position(2) = 0.99;
hText(1).String = 'ver-';
hText(1).Position(2) = 0.95;
hText(6) = copyobj(hText(1), gca);
hText(6).String = 'sal';
hText(6).Position(2) = 0.91;

% Remove legend
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend, 'Visible', 'off');
