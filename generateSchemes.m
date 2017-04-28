function generateSchemes = generateSchemes(grads,DELTAS,deltas,pGamma,...
    time180,gradMax,TE,timeConstraint,bValMin,bValMax,incFloatingPointCmp)
%inputs: list of gradients - grads
%        list of small del - deltas
%        list of big DEL - DELTAS
%        proton gamma - pGamma
%        gradMax - maximum gradient strength
%        time180 - time for 180 pulse
%        TE - echo time, which along with timeConstraint constrains DEL+del
%        timeConstraint - time for readout and period before first gradient
%        bValMin - minimum b-value
%        bValMax - maximum b-value (added 29/01/16)
%        incFloatingPointCmp - added 26/08/15
%                            - use 1 if code accounts for floating point
%                            comparisons
%                            - use 0 if floating point issues are ignored
%                            (i.e. original code)

%output: all combinations of grads, deltas and DELTAS that obey constraints

%%
%reshape inputs to ensure they are in the correct format for combining
grads = reshape(grads,1,size(grads(:),1));
DELTAS = reshape(DELTAS,1,size(DELTAS(:),1));
deltas = reshape(deltas,1,size(deltas(:),1));

%%
%first get all combinations, not all of which will be allowed
allCombinations = combvec(grads,DELTAS,deltas,pGamma)';
%assign to another variable from which certain combinations can be removed
a = allCombinations;

%%
%now get rid of unphysical combinations where delta>DELTA...
if incFloatingPointCmp==0
    a(a(:,3)>a(:,2),:)=[];
    %...and check that the indexing has worked
    if unique(a(:,3)>a(:,2))~=0
        disp(['!!! Error with scheme file filter - '...
            'unphysical combinations remain']);
    else
        %if we're here, the unphysical combinations have been removed
    end
elseif incFloatingPointCmp==1
    %27/08/15 - add additional check to account for floating point issues 
    a(a(:,3)>a(:,2) & abs(a(:,3)-a(:,2))>eps,:)=[];
    %...and check that the indexing has worked
    if unique(a(:,3)>a(:,2) & abs(a(:,3)-a(:,2))>eps)~=0
        disp(['!!! Error with scheme file filter - '...
            'unphysical combinations remain']);
    else
        %if we're here, the unphysical combinations have been removed
    end
else
    disp('!!! Error - incFloatingPointCmp must be 0 or 1')
end
%%
%now remove combinations where there is insufficient time for the 180 pulse
if incFloatingPointCmp==0
    %ie. remove if delta+time180>DELTA...
    a(a(:,3)+time180>a(:,2),:)=[];
    %...and check indexing has worked
    if unique(a(:,3)+time180>a(:,2))~=0
        disp(['!!! Error with scheme file filter - '...
            'violation of 180 pulse time constraint']);
    else
        %if we're here, we're ok
    end
elseif incFloatingPointCmp==1
    %27/08/15 - add additional check to account for floating point issues
    a(a(:,3)+time180>a(:,2) & abs(((a(:,3))+(time180))-(a(:,2)))>eps,:)=[];
    %...and check indexing has worked
    if unique(a(:,3)+time180>a(:,2) & ...
            abs(((a(:,3))+(time180))-(a(:,2)))>eps)~=0
        disp(['!!! Error with scheme file filter - '...
            'violation of 180 pulse time constraint']);
    else
        %if we're here, we're ok
    end
else
    disp('!!! Error - incFloatingPointCmp must be 0 or 1')
end
%%
%now remove combinations where delta+DELTA is too long given
if incFloatingPointCmp==0
    %ie. remove if delta+DELTA>TE-timeConstraint...
    a(a(:,3)+a(:,2)>TE-timeConstraint,:)=[];
    %...and check indexing has worked
    if unique(a(:,3)+a(:,2)>TE-timeConstraint)~=0
        disp(['!!! Error with scheme file filter - '...
            'delta+DELTA is too long for TE and timeConstraint']);
    else
        %if we're here, we're ok
    end
elseif incFloatingPointCmp==1
    %27/08/15 - add additional check to account for floating point issues
    a(a(:,3)+a(:,2)>TE-timeConstraint & ...
        abs(a(:,3)+a(:,2)-(TE-timeConstraint))>eps,:)=[];
    %...and check indexing has worked
    if unique(a(:,3)+a(:,2)>TE-timeConstraint & ...
        abs(a(:,3)+a(:,2)-(TE-timeConstraint))>eps)~=0
        disp(['!!! Error with scheme file filter - '...
            'delta+DELTA is too long for TE and timeConstraint']);
    else
        %if we're here, we're ok
    end
else
    disp('!!! Error - incFloatingPointCmp must be 0 or 1')
end

%%
%now remove combinations where gradient strength exceeds Gmax...
if incFloatingPointCmp==0
    a(a(:,1)>gradMax,:)=[];
    %...and check indexing has worked
    if unique(a(:,1)>gradMax)~=0
        disp(['!!! Error with scheme file filter - '...
            'violation of Gmax constraint']);
    else
        %if we're here, we're ok
    end
elseif incFloatingPointCmp==1
    %27/08/15 - add additional check to account for floating point issues
    a(a(:,1)>gradMax & abs(a(:,1)-gradMax)>eps,:)=[];
    %...and check indexing has worked
    if unique(a(:,1)>gradMax & abs(a(:,1)-gradMax)>eps)~=0
        disp(['!!! Error with scheme file filter - '...
            'violation of Gmax constraint']);
    else
        %if we're here, we're ok
    end
else
    disp('!!! Error - incFloatingPointCmp must be 0 or 1')
end
%%
%now remove combinations where bvalue is lower than bMin...
% 29/01/16 - also remove combinations where bvalue is greater than bMax
calc_bVal=calculate_b_value(0,a(:,1).*1e3,a(:,3).*1e3,0,a(:,2).*1e3);
if incFloatingPointCmp==0
    a(calc_bVal<bValMin | calc_bVal>bValMax,:)=[];
    %...and check indexing has worked
    if unique(calc_bVal<bValMin)~=0 | unique(calc_bVal>bValMax)~=0
        disp(['!!! Error with scheme file filter - '...
            'violation of bValMin and/or bValMax constraint']);
    else
        %if we're here, we're ok
    end

elseif incFloatingPointCmp==1    
    %27/08/15 - add additional check to account for floating point issues
    a(calc_bVal<bValMin & abs(calc_bVal-bValMin)>eps ...
      | calc_bVal>bValMax & abs(calc_bVal-bValMax)>eps,:)=[];
    %...and check indexing has worked
    if unique(calc_bVal<bValMin & abs(calc_bVal-bValMin)>eps ...
            | calc_bVal>bValMax & abs(calc_bVal-bValMax)>eps)~=0
        disp(['!!! Error with scheme file filter - '...
            'violation of bValMin and/or bValMax constraint']);
    else
        %if we're here, we're ok
    end
        
else
    disp('!!! Error - incFloatingPointCmp must be 0 or 1')
end

%%
%return a
generateSchemes=a;