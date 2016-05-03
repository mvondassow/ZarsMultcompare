%% Simulations run to test multiple comparison approach.
% These simulations were run to test whether Zar's approach to multiple
% comparisons (using Tukey's HSD) has been correctly implemented (in 
% multcompareRandom()), and gives results close to the desired alpha level
% (or CIs close to the desired nominal confidence interval) under various
% circumstances. 
%
% Ideally, the observed alpha should match the nominal alpha, and the
% fraction of time all CIs contain the 'true' difference should be 1-alpha
% (the familywise error rate should be alpha).
%
% The data generator for the simulations assumes the following model:
% Y(ijk) = F1(i) + F2(j) + F1xF2(ij) + Slope*T(ijk) + error(ijk), 
% with F1 = 0 (the null hypothesis is true for factor 1, the factor whose
% levels are compared), F2~N(0, randsd), F1xF2~N(0, interaction), and
% error~N(0,1). Setting randsd = 1 or interaction = 1 means the variance
% components associated with these terms is 1; setting slope = 1 makes 
% maximum the difference among subgroup members due to T reach 1.

clear all; close all hidden
disp(['Date run: ', date]);

%Initial values
nreps = 5000; % Number of repetitions of simulation
npergrp = 5; % Maximum sample size per subgroup (factor1xfactor2)
ntris = 4; % Number of levels of random factor (factor2)
nlvls = 3; % Number of levels of fixed factor (factor1), whose levels are
    % compared by multcompareRandom (factor 1 follows null hypothesis of no 
    % effect in simulated data).
randsd = 0; % SD of means of factor 2: (variance component)^0.5 
interaction =0; % SD of means of factor1xfactor2 interaction: 
    % (variance component)^0.5 
slope = 0; % Magnitude of effect of linear covariate term (slope = 1 would
    % give a difference of 1 between 1st and last member of each subgroup)
missing = 0; % Probability of eliminating each member of subgroup (on
    % average, max n per subgroup will be npergrp, min n per subgroup set 
    % to 1, average n per subgroup will be about (1-missing)*npergrp
alpha = 0.05; % Alpha level for comparisons (CI's set to 1-alpha).
nway = 3; % Number of 

%% Balanced 3 way ANOVA; no factors have effect. 
% no factors have effect.
testmultcompare('nreps', nreps, 'npergrp', npergrp, 'ntris', ntris, ...
    'nlvls', nlvls, 'randsd', randsd, 'interaction', interaction, ...
    'slope', slope, 'missing', missing, 'alpha', alpha, 'nway', nway);

%% Balanced 3 way ANOVA; all other factors have effects. 
% Factor2, factor1xfactor2 interaction, and linear covariate all have
% effects (factor 1 always has no effect).
testmultcompare('nreps', nreps, 'npergrp', npergrp, 'ntris', ntris, ...
    'nlvls', nlvls, 'randsd', 1, 'interaction', 1, ...
    'slope', 1, 'missing', missing, 'alpha', alpha, 'nway', nway);

%% Unbalanced 3 way ANOVA; all other factors have effects
% Factor2, factor1xfactor2 interaction, and linear covariate all have
% effects (factor 1 always has no effect).
testmultcompare('nreps', nreps, 'npergrp', npergrp, 'ntris', ntris, ...
    'nlvls', nlvls, 'randsd', 1, 'interaction', 1, ...
    'slope', 1, 'missing', 0.2, 'alpha', alpha, 'nway', nway);

%% Unbalanced 2 way ANOVA; all other factors have effects.
% Factor2 & factor1xfactor2 have effects (factor 1 always has no effect);
% covariate has no effect.
testmultcompare('nreps', nreps, 'npergrp', npergrp, 'ntris', ntris, ...
    'nlvls', nlvls, 'randsd', 1, 'interaction', 1, ...
    'slope', 0, 'missing', 0.2, 'alpha', alpha, 'nway', 2);

%% Balanced complete block without replication (2 way ANOVA). 
% Factor2 & factor1xfactor2 have effects (although interaction is included
% in error term in ANOVA analysis (factor 1 always has no effect);
% covariate has no effect.
testmultcompare('nreps', nreps, 'npergrp', 1, 'ntris', ntris, ...
    'nlvls', nlvls, 'randsd', 1, 'interaction', 1, ...
    'slope', 0, 'missing', 0, 'alpha', alpha, 'nway', 2);

%% One way ANOVA. All effects set to zero.
testmultcompare('nreps', nreps, 'npergrp', 10, 'ntris', 1, ...
    'nlvls', nlvls, 'randsd', 0, 'interaction', 0, ...
    'slope', 0, 'missing', 0, 'alpha', alpha, 'nway', 1);