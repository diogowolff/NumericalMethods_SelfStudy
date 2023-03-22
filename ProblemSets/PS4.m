clear

addpath('../AuxiliaryMethods');

%% item a)

rho = .9;
sigma = .01;

[ShockGrid, MarkovMatrix] = TauchenDiscretizer(9, 3, 0, rho, sigma);

%% item b)
beta=.96;
Value = ExperimentalHuggett(...
    [.96, 1.0001, .9, .01],  .0388, 500, 9);

%% item c)

Equilibrium = HuggettEquilibrium([.96, 1.0001, .9, .01], 500, 9);

[EqValue, UtilityMatrix, MarkovMatrix, AssetGrid] = ExperimentalHuggett(...
    [.96, 1.0001, .9, .01],  Equilibrium, 500, 9);
[~, EqInd] = PolicyHuggett(.96, EqValue, AssetGrid, 9, UtilityMatrix, MarkovMatrix);
EqDist = InvariantDistribution(EqInd, MarkovMatrix);

EarPlotGrid = linspace(-0.01, 1/beta-1, 20);
for i=1:20
    graph(i) = AssetDemandFunction(EarPlotGrid(i), Value, [.96, 1.0001, .9, .01]);
    disp(i)
end

plot(graph, EarPlotGrid);

%% Trying to replicate figure 18.6.3 from L-S

% The statement 'stdev of eps_t is .4*sqrt(1-.2^2)' from the book is
% wrong; in order to replicate, i had to set stdev of eps_t to .4

Value3 = LSHuggett([.96, 3, .2, .4], .026, 500, 7, 3);

EarPlotGrid3 = linspace(-0.01, 1/beta-1, 20);
for i=1:20
    graph3(i) = LSDemand(EarPlotGrid3(i), Value3, [.96, 3, .2, .4], 3);
    disp(i)
end

Eq3 = LSEquilibrium([.96, 3, .2, .4], 500, 7, 3);

Value6 = LSHuggett([.96, 3, .2, .4], .026, 500, 7, 6);

EarPlotGrid6 = linspace(-0.01, 1/beta-1, 20);
for i=1:20
    graph6(i) = LSDemand(EarPlotGrid6(i), Value6, [.96, 3, .2, .4], 6);
    disp(i)
end

Eq6 = LSEquilibrium([.96, 3, .2, .4], 500, 7, 6 );