clear

addpath('../AuxiliaryMethods');

%% item a)

rho = .9;
sigma = .01;

[ShockGrid, MarkovMatrix] = TauchenDiscretizer(9, 3, 0, rho, sigma);

%% item b)

[Value, UtilityMatrix, MarkovMatrix, AssetGrid] = HowardHuggett([.96, 1.0001, .9, .01], .05, 1000, 9);

%% item c)

[test, test2] = PolicyHuggett(.96, Value, AssetGrid, 9, UtilityMatrix, MarkovMatrix);
[test4, test5, test6] = InvariantDistribution(test2, MarkovMatrix);

 = HuggettEquilibrium([.96, 1.0001, .9, .01], 500, 9);

plot(test6)
for i=1:500
    test7(i) = AssetDemandFunction(i/2000, Value);
end