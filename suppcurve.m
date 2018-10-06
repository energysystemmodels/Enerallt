% This function creates a price-dependent supply curve 
% The inputs are the matrix of power supply in a region ("Supply" should be in the form of a cummulative sum),... 
% cost for each supply mode ("Cost"), and the cost range as a multiplier to base ("range", by which the supply curve is going to be adjusted)

% The output is a price-dependent supply curve ("Cost_final")

function Cost_final = suppcurve( Supply, Cost, range)

%------------------------------------------------------------------------------------
% Example:
% You can convert this function to script and see how it works with input data below
% Supply=1.0e+04 *[0.0010 0.0357 0.2642 0.2651 0.5403 0.7088 0.9139 0.9211 1.2003 1.3614 1.5319];
% Cost=[ 2       4           10   15 30   35    40  45   55   65  75];
% range=0.0;
%-------------------------------------------------------------------------------------

Supply=fix(Supply);
Cost_final = zeros(1,max(Supply));

if length(Cost)==1                           % If only one production mode (one cost)
  for i=1:Supply(1)
  Cost_final(i) = Cost(1)+(i-0)*(range*(Cost(1)-0))/(Supply(1)-0);
  end
else

% If more than one production mode

for i=1:Supply(1)
  Cost_final(i) = Cost(1)+(i-0)*(range*(Cost(2)-Cost(1)))/(Supply(1)-0);
end
for j=1:length(Supply)-2
     for i=Supply(j)+1:Supply(j+1)
      Cost_final(i) =Cost(j+1)+(i-Supply(j))*(range*(Cost(j+2)-Cost(j+1)))/(Supply(j+1)-Supply(j));
    end
end
for i=Supply(end-1)+1:Supply(end)
Cost_final(i) = Cost(end)+(i-Supply(end-1))*(range*(Cost(end)-Cost(end-1)))/(Supply(end)-Supply(end-1));
end
end

