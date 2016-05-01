function [cimatch, pmatch] = testmultcompare(varargin)
%Simulation to test multiple comparison procedure using a 3 way ANOVA with
%possible unbalanced design
% Inputs as name value pairs:
%   nreps : integer, number of repetitions of simulation
%   ngrp : integer, (maximum) number of individuals per subgroup
%   missing : fraction

% Default settings
nreps = 10;
npergrp = 5;
ntris = 4;
nlvls = 3;
factor2 = 0;
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
        case 'factor2'
            factor2 = varargin{i+1};
        case 'covariate'
            covarate = varargin{i+1};
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
    simtr=repmat(repmat(rand(1,1)*randn(ntris,1), nlvls, 1), 1, npergrp);
    %Simulate trialxtreatment interactions, with magnitude factor2
    simx=repmat(factor2*randn(ntris*nlvls,1), 1, npergrp);
    
    %simulate time effect with magnitude 'covariate'
    tlist=repmat(0.5:npergrp-0.5, nlvls*ntris, 1)-npergrp/2;
    teffect=tlist*covariate;
    
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
    if size(grpstats(simdat2, {treatlist2, trilist2}), 1)==12      
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
        
        
    end
end

disp(['alpha: ', num2str(alpha)]);
disp(['type 1 error rate: ', ...
    num2str(sum(sigcomp)/numel(sigcomp))]);
disp(['proportion fail to contain true dif: ',...
    num2str(sum(outofci)/numel(outofci))]);
disp(['proportion ANOVA found factor 1 significant: ', ...
    num2str(sum(siganova)/numel(siganova))]);

end

