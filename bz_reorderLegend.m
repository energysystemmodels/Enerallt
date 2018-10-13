function out = bz_reorderLegend(order,h)
% Reorders the legend entries (without modifying the actual order the curves are stacked on the axes)
% order  - an array specifying the order of legend entries
% h - axis handle (optional, defaults to current axis)
    if nargin < 2
        h = get(gca,'Children');
    end
     
    hC = findobj(h);
    
    if nargin == 0 || isempty(order)
        % output/display the internal order
        lbl = get(hC,{'DisplayName'});        
        if nargout == 1 
            out = lbl;
        else
            for k = 1:length(lbl)
               display([sprintf('[%*d] ',fix(log10(length(lbl)))+1,k) lbl{k}]);
            end
        end
        
    else
       % reorder the legend entries 
       legend(hC(order)); 
    end
end