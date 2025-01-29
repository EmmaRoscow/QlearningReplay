
colour_scheme;
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
ksamples = [0 1 3 5 10 15 20 30 40 50 75 100];

% Colourmap: black to blue to purple
colour_anchors = [0, 0, 0;...
                         1, 0.1, 0.1;...
                         colourscheme.RPEProb];

% Black to red transition
cmap1 = [linspace(colour_anchors(1, 1), colour_anchors(2, 1), length(ksamples)*10/2)',...
    linspace(colour_anchors(1, 2), colour_anchors(2, 2), length(ksamples)*10/2)',...
    linspace(colour_anchors(1, 3), colour_anchors(2, 3), length(ksamples)*10/2)'];

% Red to orange transition
cmap2 = [linspace(colour_anchors(2, 1), colour_anchors(3, 1), length(ksamples)*10/2)',...
    linspace(colour_anchors(2, 2), colour_anchors(3, 2), length(ksamples)*10/2)',...
    linspace(colour_anchors(2, 3), colour_anchors(3, 3), length(ksamples)*10/2)'];

% Join
cmap = [cmap1; cmap2];

% Plot
[h] = plot_reliability_error_over_sessions_x_replays(rat_names, ksamples, 'RPEProb', cmap);

% Adjust axis
yyaxis left;
ylim([0.42 2.9])


