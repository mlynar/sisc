function cm = com(x, positiveOnly)
%index of center of mass
%x - column vector
%positiveOnly - if true returns pos of max(x) if false returns pos of
%max(abs(x))
    %trs = 1e-3;
    %x = x(abs(x) > trs);
    
    %L = length(x);
    %cm = floor(((1:L) * abs(x)));
    
    if nargin < 2
        positiveOnly = true;
    end
    
    if positiveOnly
        cm = find(x == max(x));
    else
        cm = find(abs(x) == max(abs(x)));
    end
