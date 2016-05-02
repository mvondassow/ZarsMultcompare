function [simresults] = testmultcompare(varargin)
%Simulation to test multiple comparison procedure using a 3 way ANOVA with
%possible unbalanced design (using name-value pair 'missing'). First
%factor is a fixed factor with no effect (follows null hyp), factor2 (or
%trial=tri) is a random factor, and covariate is a linear covariate.
%factor1 and factor2 can have an interaction (set by 'interaction' name
%value pair), but the covariate does not. e.g. 
%Y(ijk) = F2(j) + F1xF2(ij) + T(ijk) +error(ijk) 
% Inputs as name-value pairs:
%   nreps : integer, number of repetitions of simulation
%   npergrp : integer, maximum number of individuals per subgroup
%   missing : float, fraction of original set to delete (to give unbalanced
%       design). 0<missing<1
%   tris : integer, number of levels of random factor.
%   nlvls : integer, number of levels of fixerd factor.
%   factor2 : magnitude of random factor effect (variance^0.5)
%   interaction : magnitude of interaction between fixed and random factors
%       (approximately variance^0.5)
%   covariate : magnitude of linear covariate, scaled so that a value of 1
%       means a difference of 1 between the first and last member of each
%       subgroup (counting missing entries).
%   NOTE: Error term is fixed as normal distribution with mean=0, SD=1.
%   alpha : alpha level (nominal error rate) for confidence intervals and 
%       pairwise tests
%
% Returns
%   cimatch : observed error rate (excluding true value of 0 for difference
%       among levels) for any of the confidence intervals for differences.
%   pmatch : observed type 1 familywise error rate for pairwise tests. 
%   cimatch and pmatch should approximate alpha.

% Default settings
nreps = 10;
npergrp = 5;
ntris = 4;
nlvls = 3;
randfactor = 0;
interaction = 0;
covariate = 0;
missing = 0;
alpha = 0.05;

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
        case 'randfactor'
            randfactor = varargin{i+1};
        case 'interaction'
            interaction = varargin{i+1};
        case 'covariate'
            covariate = varargin{i+1};
        otherwise
            error('Invalid name-value pair')
    end
end

% initials arrays to keep instances identified a significant differences 
% or failure to contain 0 in CI
sigcomp = nan(nreps, 1);
outofci = nan(nreps, 1);
siganova = nan(nreps, 1); % Check that ANOVA is working as expected.

for k=1:nreps
%Generate random data:
    %Generate labels for simulated treatments and trials
    treatlist = repmat(...
        reshape(repmat([1:nlvls], ntris, 1), nlvls*ntris, 1), 1, npergrp);
    trilist = repmat(repmat([1:ntris]', nlvls,1),1, npergrp);
    %Simulate trial effect, with randomized magnitude;
    simtr=repmat(repmat(randfactor*randn(ntris,1), nlvls, 1), 1, npergrp);
    %Simulate trialxtreatment interactions, with magnitude randfactor
    simx=repmat(interaction*randn(ntris*nlvls,1), 1, npergrp);
    
    %simulate covariate (e.g. time) effect with magnitude 'covariate'
    tlist=repmat(0.5:npergrp-0.5, nlvls*ntris, 1)-npergrp/2;
    teffect=tlist*covariate/(npergrp-1);
    
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
    if size(grpstats(simdat2, {treatlist2, trilist2}), 1)==nlvls*ntris      
        %run anovan on simulated data to gets stats table.
        [Panova,~,simstats]=anovan(...
            simdat2, {treatlist2, trilist2, tlist2}, ...
            'random', 2, 'model', [1,0,0;0,1,0;1,1,0;0,0,1], ...
            'continuous', 3, 'display', 'off');
     
        %Run multiple comparison
        [Pmat, CImat] = multcompareRandom(...
            simdat2, treatlist2, simstats, alpha, false);
        
        % Test if any comparisons are found to be significantly different.
        sigcomp(k) = any(Pmat<alpha);
        % Test if any CIs fail to contain 0.
        outofci(k) = any([(CImat(:,1)>0); (CImat(:,3)<0)]);
        siganova(k) = (Panova(1)<alpha);
                
    else
        %Enter elimnatable value for repetitions in which some
        %treatment*trial subgroup was eliminated.
        sigcomp(k) = -1;
        outofci(k) = -1;
        siganova(k) = -1;
    end
end
successes = (sigcomp~=-1);
sigcomp=sigcomp(successes);
outofci=outofci(successes);
siganova=siganova(successes);

disp(['Nominal type 1 error rate (alpha): ', num2str(alpha)]);
disp(['Observed type 1 error rate: ', ...
    num2str(sum(sigcomp)/numel(sigcomp))]);
disp(['Proportion of times that any CIs fail to contain true dif (0): ',...
    num2str(sum(outofci)/numel(outofci))]);
disp(['Observed type 1 error rate for ANOVA: ', ...
    num2str(sum(siganova)/numel(siganova))]);
disp(['Number of completed iterations: ', num2str(numel(sigcomp))]);

simresults.sigcomp=sigcomp;
simresults.outofci=outofci;
simresults.siganova=siganova;
end

