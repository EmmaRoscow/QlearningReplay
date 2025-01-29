
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
policies = {'random', 'expectedReward', 'RPE', 'RPEProb'};
shuffled = true;

% Plot
[h, t] = plot_reliability_error_over_replays(rat_names, policies, shuffled);

% Adjust legend
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.Position(1) = hLegend.Position(1)-0.01;
hLegend.Position(2) = hLegend.Position(2)+0.03;