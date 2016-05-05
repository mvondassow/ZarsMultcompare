function [multps, cimat] = multcompareRandom(myvar, ...
    factor1, mystats, varargin)
% Do multiple comparisons and CI procedure from Zar (Tukey's method).
% 
% Assumes either a 1 way ANOVA or a multiway ANOVA, for multiway ANOVAs
% (2-way, 3-way, etc) first factor must be a fixed factor whose levels will
% be compared, and second factor must be a categorical random factor.
%
% Inputs (positional):
%   myvar : numeric 1D array of dependent variable
%   factor1 : 1D array listing treatment levels of factor 1 for each datum
%   mystats : stats structure produced by ANOVAN
%   alpha : (optional) float<1. Alpha value for comparison (default = 0.05)
%   verbose : (optional) 0, true or 0, false. Whether to display a table
%       with pairwise comparisons, CIs, and P values.
%
% Outputs: 
%   multps : array with P values for pairwise comparisons.
%   citmat : array with CIs for differences among treatment levels for
%       factor 1 (first column: lower bound; middle column: mean; third
%       column: upper bound)

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
    
    % Check whether it is one way or multi-way anova, and identify msdenom 
    % and dfdenom accordingly (based on how MATLAB generates stats
    % structure).
    if isequal(mystats.terms, 1)
        msdenom = mystats.mse;
        dfdenom = mystats.dfe;
    else
        % if it is multiway ANOVA check that it has first factor fixed, 
        % and second factor random.
        if ~strcmp(mystats.varnames(2),mystats.rtnames(1))
            error('First factor must be fixed; second factor must be random');
        end
        dfdenom = mystats.dfdenom(1);
        msdenom = mystats.msdenom(1);
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
    
    % Calculate SE following Zar. 
    se = (msdenom.*...
            (1./grpns(mypairs(:,1))+1./grpns(mypairs(:,2)))/2).^0.5;
        
    %Calculate p-values for differences in group means.
    multps=nan(nl1,1);
    for l=1:mystats.nlevels(1)
       multps(l)=1-stdrcdf(...
             abs(mt(mypairs(l,1))-mt(mypairs(l,2)))/se(l),...
             round(dfdenom), nl1);
    end
    
    %Calculate 1-alpha ci's for differences among means.
    cidif=stdrinv(1-alpha,round(dfdenom),nl1).*se;
    %Estimate of mean difference:
    meandifs=mt(mypairs(:,1))-mt(mypairs(:,2));
    
    %Calculate confidence intervals for differences (column 1: lower
    %bounds, column 2: estimated mean differences, column 3: upper bounds)
    cimat = [meandifs-cidif, meandifs, ...
            meandifs+cidif];
    
    if verbose == 1
        %create table for CI's and multiple comparisons
        disp([[{'treatment 1'}, {'treatment 2'}, ...
            {['LB ', num2str(100*(1-alpha)), '% CI']},...
            {'mean difference'}, ...
            {['UB ', num2str(100*(1-alpha)), '% CI']}, {'P, pairwise'}];...
            [treatname(mypairs(:,1)), treatname(mypairs(:,2)),...
            num2cell(cimat), num2cell(multps)]]);

        %Create table of estimated treatment means (estimated from coeffs
        %calculated by ANOVAN).
        disp([[{'Treatment'}, {'Mean'}];[treatname(:), num2cell(mt)]]);
    end

end

