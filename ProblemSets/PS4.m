clear

addpath('../AuxiliaryMethods');

%% item a)

rho = .9;
sigma = .01;

[ShockGrid, MarkovMatrix] = TauchenDiscretizer(9, 3, 0, rho, sigma);

%% item b)

[Value, ~, ~, AssetGrid] = ExperimentalHuggett(...
    [.96, 1.0001, .9, .01],  .03, 500, 9);

set(gcf,'DefaultLineLineWidth', 1.5);
plot(AssetGrid, Value);
title('Value function with r = 0.03');
xlabel('asset level, $a$','Interpreter','latex');
ylabel('value function, $v(a,\bar{y})$','Interpreter','latex');
ha = legend('$y=y_1$', '$y=y_2$', '$y=y_3$',...
            '$y=y_4$', '$y=y_5$', '$y=y_6$',...
            '$y=y_7$', '$y=y_8$', '$y=y_9$',...
            'Location','SouthEast');
set(ha,'Interpreter','latex');

%% item c)

Equilibrium = HuggettEquilibrium([.96, 1.0001, .9, .01], 500, 9);

[EqValue, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett(...
    [.96, 1.0001, .9, .01],  Equilibrium, 500, 9);
[~, EqInd] = PolicyHuggett(.96, EqValue, AssetGrid, 9, UtilityMatrix, MarkovMatrix);
EqDist = InvariantDistribution(EqInd, MarkovMatrix);

[X,Y] = meshgrid(AssetGrid,exp(ShockGrid));
surf(X,Y, EqDist', EqValue');
hcb = colorbar;
xlabel('asset level, $a$','Interpreter','latex');
ylabel('endowment, $y$','Interpreter','latex');
zlabel('measure $\mu$','Interpreter','latex');
title('Invariant distribution at equilibrium');
colorTitleHandle = get(hcb,'Title');
titleString = 'Value function';
set(colorTitleHandle ,'String',titleString);


EarPlotGrid = linspace(-0.01, 1/beta-1, 20);
for i=1:20
    graph(i) = AssetDemandFunction(EarPlotGrid(i), EqValue, [.96, 1.0001, .9, .01]);
    disp(i)
end

plot(graph, EarPlotGrid);
xlabel('Expected demand \bf{E}[a(r)]');
ylabel('interest rate, r');
title('Aggregate demand of assets');
xline(0);
yline(0);


%% item d)

Equilibriumd = HuggettEquilibrium([.96, 1.0001, .97, .01], 200, 9);

%% item e)

Equilibriume = HuggettEquilibrium([.96, 5, .9, .01], 200, 9);

%% item f)

Equilibriumf = HuggettEquilibrium([.96, 1.0001, .9, .05], 500, 9);

EarPlotGridf = linspace(-0.01, 1/beta-1, 20);
for i=1:20
    graphf(i) = AssetDemandFunction(EarPlotGridf(i), EqValue, [.96, 1.0001, .9, .05]);
    disp(i)
end

plot(graph, EarPlotGrid);
hold on;
plot(graphf, EarPlotGridf);
legend('Original model', 'Model with higher variance',...
        'Location', 'SouthEast');

%% Trying to replicate figure 18.6.3 from L-S

% The statement 'stdev of eps_t is .4*sqrt(1-.2^2)' from the book is
% wrong; in order to replicate, i had to set stdev of eps_t to .4

Value3 = LSHuggett([.96, 3, .2, .4], .026, 500, 7, 3);

EarPlotGrid3 = linspace(-0.01, .041, 20);
for i=1:20
    graph3(i) = LSDemand(EarPlotGrid3(i), Value3, [.96, 3, .2, .4], 3);
    disp(i)
end

%Eq3 = LSEquilibrium([.96, 3, .2, .4], 500, 7, 3);

Value6 = LSHuggett([.96, 3, .2, .4], .026, 500, 7, 6);

EarPlotGrid6 = linspace(-0.01, .041, 20);
for i=1:20
    graph6(i) = LSDemand(EarPlotGrid6(i), Value6, [.96, 3, .2, .4], 6);
    disp(i)
end

plot(graph3, EarPlotGrid3);
hold on;
plot(graph6, EarPlotGrid6);
xline(0);
yline(0);
title('Graph 18.6.3')

%Eq6 = LSEquilibrium([.96, 3, .2, .4], 500, 7, 6 );

%% Figure 18.7.2 invariant distribution

[Eq3] = LSEquilibrium([.96, 3, .2, .4], 500, 7, 3);

[EqValue3, UtilityMatrix3, MarkovMatrix3, AssetGrid3] = LSHuggett(...
    [.96, 3, .2, .4], 0, 100, 7, 3);
[~, EqInd3] = PolicyHuggett(.96, EqValue3, AssetGrid3, 7, UtilityMatrix3, MarkovMatrix3);
EqDist3 = InvariantDistribution(EqInd3, MarkovMatrix3);

plot(AssetGrid3, sum(EqDist3, 2))
