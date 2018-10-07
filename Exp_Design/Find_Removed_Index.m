function [Removed_States,Removed_Index] = Find_Removed_Index(Num_var)
% This function finds the index of the functions who are not included in
% the search spave of Num_var input variables because they do not include
% all Num_var input variables in their logics


load('Target_Output_Full.mat');

Num_State = 2^Num_var;

Input_States = de2bi(Num_State-1:-1:0,'left-msb');

%Remove all 0 and 1, and A,~A, ....
Removed_States = [zeros(Num_State,1) ones(Num_State,1)];
for i=1:Num_var
    Removed_States = [Removed_States Input_States(:,i) ~Input_States(:,i)];
end
input_Index = Num_State- bi2de(Input_States,'left-msb');

switch Num_var
    
    case 3
        
        Indexes = [];
        for ii=1:3
            for jj=1:3
                if jj==ii
                    continue;
                end
                for kk=1:3
                    if(kk==ii || kk==jj)
                        continue;
                    end
                    Input_Order = [ii jj kk];
                    [~,i_ind]=sort([ii jj kk]);
                    Reorder_Index = Num_State- bi2de(Input_States(:,i_ind),'left-msb');
                    [~, ind]=sort(Reorder_Index);
                    Indexes = [Indexes ind];
                end
            end
        end
        
        
        %Find target outputs of (A,C), (A,C) and (B,C)
        Removed_States = Calculate_Removed_States(Num_var,Target_Output,Removed_States,Indexes);
        
    case 4
        Indexes = [];
        for ii=1:4
            for jj=1:4
                if jj==ii continue;  end
                for kk=1:4
                    if(kk==ii || kk==jj)  continue;  end
                    for mm=1:4
                        if(mm==ii || mm==jj || mm==kk)  continue;  end
                        
                        Input_Order = [ii jj kk mm];
                        [~,i_ind]=sort([ii jj kk mm]);
                        Reorder_Index = Num_State- bi2de(Input_States(:,i_ind),'left-msb');
                        [~, ind]=sort(Reorder_Index);
                        Indexes = [Indexes ind];
                    end
                end
            end
        end
        
        
        %Find target outputs of (A,C), (A,C) and (B,C)
        Removed_States = Calculate_Removed_States(Num_var,Target_Output,Removed_States,Indexes);
        
        
    case 5
        Indexes = [];
        for ii=1:5
            for jj=1:5
                if jj==ii continue;  end
                for kk=1:5
                    if(kk==ii || kk==jj)  continue;  end
                    for mm=1:5
                        if(mm==ii || mm==jj || mm==kk)  continue;  end
                        for nn=1:5
                            if(nn==ii || nn==jj || nn==kk || nn==mm)  continue;  end
                            
                            Input_Order = [ii jj kk mm nn];
                            [~,i_ind]=sort([ii jj kk mm nn]);
                            Reorder_Index = Num_State- bi2de(Input_States(:,i_ind),'left-msb');
                            [~, ind]=sort(Reorder_Index);
                            Indexes = [Indexes ind];
                        end
                    end
                end
            end
        end
        
        
        %Find target outputs of (A,C), (A,C) and (B,C)
        Removed_States = Calculate_Removed_States(Num_var,Target_Output,Removed_States,Indexes);
end

Removed_States = unique(Removed_States','rows')';

%The code is for reverse direction of input states (all one at top and all
%0 at the buttom. We do flipud to move all 0 to top
Removed_States = flipud(Removed_States);

%Find the function index from output states
Removed_Index = bi2de(Removed_States','left-msb')';

%Sort the indexes
[Removed_Index,indx] = sort(Removed_Index);
Removed_States = Removed_States(:,indx);

end