% This function creates a price-dependent supply curve 
% The inputs are the matrix of power supply (it should be in the form of
% cummulative sum), cost for each supply mode, and the coefficient by which
% the supply curve is going to be adjusted

function [Sup_dis,Cost_dis] = suppcurve_dis( Sup, Cost, coff, n)

    Sup_dis=zeros(1,length(Cost)*n);
    Cost_dis=zeros(1,length(Cost)*n);


    if length(Cost)==1
      for j=1:n
          Cost_dis(j)=Cost(1)+j*(coff/n)*(Cost(1));

      end
      Sup_dis(1:n)=linspace(0,Sup(1),n);

    else

    % Discrete supply curve
        for j=1:n
              Cost_dis(j)=Cost(1)+j*(coff/n)*(Cost(1));
        end
        Sup_dis(1:n)=linspace(0,Sup(1),n);

        for i=2:length(Cost)-1
            for j=1:n
                Cost_dis((i-1)*n+j)=Cost(i)+j*(coff/n)*(Cost(i+1)-Cost(i));

            end
            Sup_dis((i-1)*n+1:i*n)=linspace(Sup(i-1),Sup(i),n);
        end

        for j=1:n
                Cost_dis((length(Cost)-1)*n+j)=Cost(end)+j*(coff/n)*(Cost(end)-Cost(end-1));

        end
        Sup_dis((length(Cost)-1)*n+1:length(Cost)*n)=linspace(Sup(end-1),Sup(end),n);

    end