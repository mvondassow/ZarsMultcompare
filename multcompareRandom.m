function [multps, cimat] = multcompareRandom(myvar, ...
    factor1, mystats, varargin)
% Do multiple comparisons and CI procedure from Zar (Tukey's method).
% Inputs:
%   myvar : numeric 1D array of dependent variable
%   factor1 : 1D array listing treatment levels of factor 1 for each datum
%   mystats : stats structure produced by ANOVAN
%   alpha : (optional) float<1, alpha value for comparison (default = 0.05)
% Outputs: 
%   citable : cell array with one header row of cells, one column naming
%       treatment pairs (strings; column 1), then lower bound of CI
%       (numeric; column 2), mean difference (numeric; column 3), upper
%       bound of CI (numeric; column 3), P value(numeric; column 4).
%   meantable : cell array with treatment manes (string; column 1) and
%       estimated means (numeric; column 2)

    % Check that arguments are valid
    if nargin > 3
        alpha = varargin{1};
        verbose = 0;
        if (nargin == 5)
            if (varargin{2}==0||varargin{2}==false)
                verbose = 0;
            elseif (varargin{2}==1||varargin{2}==true)
                verbose = 1;
            else
                error('Unexpected inputs for verbose (0/1)');
            end
        elseif nargin > 5
            error('Unexpected number of inputs');   
        end
    else
        verbose = 0;
        alpha = 0.05;
    end
    
    % Check that is a random model
    if ~strcmp(mystats.varnames(2),mystats.rtnames(1))
        error('First factor must be fixed; second factor must be random');
    end

    % Names and sample size per subgroup (at level of treatment)
    [grpns, grpnames] = grpstats(myvar, factor1, {'numel', 'gname'});
    
    %Get number of levels of factor 1
    nl1=mystats.nlevels(1);
    
    % Calculate means for treatments (mt) based on coeffs in stats
    % structure from ANOVAN
    mt=mystats.coeffs(1)+mystats.coeffs(2:1+nl1);
    
    %Get treatment names
    treatname = strrep(mystats.coeffnames(2:1+nl1),...
        [mystats.varnames{1},'='],'');
    %Verify treatment names
    if min(strcmp(treatname, grpnames))==0
        error('Problem identifying treatment level names')
    end
    
    %mypairs contains a list of pairs of treatments for the first factor.
    mypairs = nchoosek([1:nl1], 2);
    
    % calculate SE. Following Zar, I based SE on the number of data points
    % per treatment and ms-error. (but seems like more sensible to use 
    % denominator ms (mystats.msdenom(1)) and n = df denominator)
    se = (mystats.msdenom(1).*...
       (1./grpns(mypairs(:,1))+1./grpns(mypairs(:,2)))/2).^0.5;
    % Was way too conservative with Zar's method.
    % Try the following:
    %se = (mystats.msdenom(1)/mystats.nlevels(2))^0.5;
    % Even more conservative

    %Calculate p-values for differences in group means.
    multps=nan(nl1,1);
    for l=1:mystats.nlevels(1)
       multps(l)=1-stdrcdf(...
             abs(mt(mypairs(l,1))-mt(mypairs(l,2)))/se(l),...
             round(mystats.dfdenom(1)), nl1);
    end
    
    %Calculate 95% ci's for differences among means.
    cidif=stdrinv(1-alpha,round(mystats.dfdenom(1)),nl1).*se;
    %Estimate of mean difference:
    meandifs=mt(mypairs(:,1))-mt(mypairs(:,2));
    
    %Calculate confidence intervals for differences (column 1: lower
    %bounds, column 2: estimated mean differences, column 3: upper bounds)
    cimat = [meandifs-cidif, meandifs, ...
            meandifs+cidif];
    
    if verbose == 1
        %create table for CI's and multiple comparisons
        disp([[{'treatment 1'}, {'treatment 2'}, {'LB 95% CI'},...
            {'mean difference'}, {'UB 95% CI'}, {'P, pairwise'}];...
            [treatname(mypairs(:,1)), treatname(mypairs(:,2)),...
            num2cell(cimat), num2cell(multps)]]);

        %Create table of estimated treatment means (estimated from coeffs
        %calculated by ANOVAN).
        disp([[{'Treatment'}, {'Mean'}];[treatname(:), num2cell(mt)]]);
    end

end

