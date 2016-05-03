function [simresults] = testmultcompare(varargin)
% Simulation to test multiple comparison procedure following ANOVA with
% possible unbalanced design (using name-value pair 'missing').
% 
% The simulated data generating process is as follows: 
% First factor is a fixed factor with no effect (it follows null hyp),
% factor2 (or trial=tri) is a random factor, and third factor is a linear
% covariate. Factors 1 and 2 can have an interaction (set by 'interaction'
% name value pair), but the covariate does not. e.g. The data generating
% model is: 
% Y(ijk) = F1(i) + F2(j) + F1xF2(ij) + Slope*T(ijk) + error(ijk),
% with F1 = 0, F2~N(0, randsd), F1xF2~N(0, interaction), and error~N(0,1)
% Data is always generated according to this model, but the this function
% can be used to test lower order models.
%
% Inputs are name-value pairs:
%   nreps : integer, number of repetitions of simulation
%   npergrp : integer, maximum number of individuals per subgroup. If
%       npergrp == 1, then cannot include covariate, or missing values; if
%       nway == 1, npergrp is still the size of the subgroups (F1xF2), even
%       though F2 is not included in the ANOVA (it is included in the data
%       generation).
%   missing : float, fraction of original set to delete (to give unbalanced
%       design). 0<missing<1
%   tris : integer, number of levels of random factor.
%   nlvls : integer, number of levels of fixerd factor.
%   randsd : magnitude of random factor effect (variance^0.5)
%   interaction : magnitude of interaction between fixed and random factors
%       (variance^0.5)
%   slope : magnitude of linear covariate (factor 3), scaled so that a
%       value of 1 gives a difference of 1 between the first and last 
%       member of each subgroup (counting missing entries).
%   NOTE: Error term is fixed as normal distribution with mean=0, SD=1.
%   alpha : alpha level (nominal error rate) for confidence intervals and 
%       pairwise tests
%   nway : integer (1, 2, or 3 only) for 1, 2, or 3 way ANOVA.
% Displays results using disp().
%
% Returns the following variables:
%   cimatch : observed error rate (excluding true value of 0 for difference
%       among levels) for any of the confidence intervals for differences.
%   pmatch : observed type 1 familywise error rate for pairwise tests. 
%   cimatch and pmatch should approximate alpha.

% Default settings
nreps = 10;
npergrp = 5;
ntris = 4;
nlvls = 3;
randsd = 0;
interaction = 0;
slope = 0;
missing = 0;
alpha = 0.05;
nway = 3;

% Switch trap parses the varargin inputs
len = length(varargin);
% check "len" for even number
if mod(len,2) > 0
    error('Wrong arguments: must be name-value pairs.');
end
for i = 1:2:len
    switch lower(varargin{i})
        case 'nreps'
            nreps=varargin{i+1};
        case 'npergrp'
            npergrp=varargin{i+1};
        case 'missing'
            missing=varargin{i+1};
        case 'alpha'
            alpha=varargin{i+1};
        case 'ntris'
            ntris = varargin{i+1};
        case 'nlvls'
            nlvls = varargin{i+1};
        case 'randsd'
            randsd = varargin{i+1};
        case 'interaction'
            interaction = varargin{i+1};
        case 'slope'
            slope = varargin{i+1};
        case 'nway'
            nway = varargin{i+1};
        otherwise
            error(['Invalid name-value pair: ', varargin{i}])
    end
end

% Generate variables to define ANOVA model
if isequal(nway, 1)
    varlist = @(x,y,z) {x};
    randfactor = [];
    covariate = [];
    mymodel = 1;
elseif isequal(nway, 2)
    varlist = @(x,y,z) {x,y};
    randfactor = 2;
    covariate = [];
    mymodel = [1,0;0,1;1,1];
    if ntris == 1
        error('For 2 or 3-way ANOVA, factor 2 #levels (ntris) must be >1)')
    end
elseif isequal(nway, 3)
    varlist = @(x,y,z) {x,y,z};
    randfactor = 2;
    covariate = 3;
    mymodel = [1,0,0;0,1,0;1,1,0;0,0,1];
    if npergrp == 1
        error('Including covariate term (nway=3) requires npergrp>1');
    end
    if ntris == 1
        error('For 2 or 3-way ANOVA, factor 2 #levels (ntris) must be >1)')
    end
else
    error('Option nway must be [1,2,3]');
end

% initials arrays to keep instances identified a significant differences 
% or failure to contain 0 in CI
sigcomp = nan(nreps, 1);
outofci = nan(nreps, 1);
siganova = nan(nreps, 1); % Check that ANOVA is working as expected.
minns = nan(nreps, 1);
maxns = nan(nreps, 1);
medianns = nan(nreps, 1);

for k=1:nreps
%Generate random data:
    %Generate labels for simulated treatments and trials
    treatlist = repmat(...
        reshape(repmat([1:nlvls], ntris, 1), nlvls*ntris, 1), 1, npergrp);
    trilist = repmat(repmat([1:ntris]', nlvls,1),1, npergrp);
    %Simulate trial effect, with randomized magnitude;
    simtr=repmat(repmat(randsd*randn(ntris,1), nlvls, 1), 1, npergrp);
    %Simulate trialxtreatment interactions, with magnitude randsd
    simx=repmat(interaction*randn(ntris*nlvls,1), 1, npergrp);
    
    %simulate covariate (e.g. time) effect with magnitude 'slope'
    tlist=repmat(0.5:npergrp-0.5, nlvls*ntris, 1)-npergrp/2;
    if npergrp>1
        teffect=tlist*slope/(npergrp-1);
    else
        teffect=tlist*0;
    end
    
    %Add random error term to each
    erlist=randn(nlvls*ntris,npergrp);
    
    %combined effects
    simdat=erlist(:)+teffect(:)+simx(:)+simtr(:);
        
    %Make gaps to simulate missing data.
    mylist=(rand(1, npergrp*nlvls*ntris)>missing);
    treatlist2=treatlist(mylist);
    trilist2=trilist(mylist);
    tlist2=tlist(mylist);
    simdat2=simdat(mylist);
    
    %Make sure I haven't eliminated any treatment/trial combinations
    nlist = grpstats(simdat2, {treatlist2, trilist2}, 'numel');
    if size(nlist, 1)==nlvls*ntris      
        %run anovan on simulated data to gets stats table.
        [Panova,simtable,simstats]=anovan(...
            simdat2, varlist(treatlist2, trilist2, tlist2), ...
            'random', randfactor, 'model', mymodel, ...
            'continuous', covariate, 'display', 'off');
     
        %Run multiple comparison
        [Pmat, CImat] = multcompareRandom(...
            simdat2, treatlist2, simstats, alpha, false);
        
        % Test if any comparisons are found to be significantly different.
        sigcomp(k) = any(Pmat<alpha);
        % Test if any CIs fail to contain 0.
        outofci(k) = any([(CImat(:,1)>0); (CImat(:,3)<0)]);
        % Test whether ANOVA identified factor 1 as significant.
        % Store information about sub-group sizes.
        siganova(k) = (Panova(1)<alpha);
        minns(k) = min(nlist);
        medianns(k) = median(nlist);
        maxns(k) = max(nlist);
                
    else
        %Enter elimnatable value for repetitions in which some
        %treatment*trial subgroup was eliminated.
        sigcomp(k) = -1;
        outofci(k) = -1;
        siganova(k) = -1;
        minns(k) = -1;
        medianns(k) = -1;
        maxns(k) = -1;
    end
end

%Eliminate entries where ANOVA wasn't run due to missing subgroup
successes = (sigcomp~=-1);
sigcomp=sigcomp(successes);
outofci=outofci(successes);
siganova=siganova(successes);

%Display information about the anova model.
disp('ANOVA MODEL');
if nway > 1
    disp([simtable(:,1), simtable(:,8)]);
    disp(['Continuous variable:',...
        simstats.varnames(simstats.continuous==1)]);
else
    disp(simtable(:,1));
end
% Display information about how unbalanced the ANOVA is.
disp(['median min subgroup n: ' num2str(median(minns(successes)))]);
disp(['median median subgroup n: ' num2str(median(medianns(successes)))]);
disp(['median max subgroup n: ' num2str(median(maxns(successes)))]);
disp(' '); disp(['Nominal type 1 error rate (alpha): ', num2str(alpha)]);

disp(' ');disp('Simulation results:');
disp(['Observed type 1 error rate: ', ...
    num2str(sum(sigcomp)/numel(sigcomp))]);
disp(['Fraction of time all CIs contain true difference: ',...
    num2str(1-sum(outofci)/numel(outofci))]);
disp(['Observed type 1 error rate for ANOVA: ', ...
    num2str(sum(siganova)/numel(siganova))]);
disp(['Number of completed iterations: ', num2str(numel(sigcomp))]);
disp(' ')

simresults.sigcomp=sigcomp;
simresults.outofci=outofci;
simresults.siganova=siganova;
end

