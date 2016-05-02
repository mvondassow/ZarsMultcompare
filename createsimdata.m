npergrp = 5;
ntris = 4;
nlvls = 3;
factor2 = 0;
interaction = 0;

covariate = 0;
missing = 0;

%Generate random data:
    %Generate labels for simulated treatments and trials
    treatlist = repmat(...
        reshape(repmat([1:nlvls], ntris, 1), nlvls*ntris, 1), 1, npergrp);
    trilist = repmat(repmat([1:ntris]', nlvls,1),1, npergrp);
    %Simulate trial effect, with randomized magnitude;
    simtr=repmat(repmat(factor2*randn(ntris,1), nlvls, 1), 1, npergrp);
    %Simulate trialxtreatment interactions, with magnitude factor2
    simx=repmat(interaction*randn(ntris*nlvls,1), 1, npergrp);
    
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