function [Gate_Cost_Table,Pr_Input_State_Mtx,Next_Initial_State,Next_Dest_State,Cost] = Find_Minimum_Cost_Transition(Pr_Input_State_Mtx,W_D,W_O,SPSD_Observed_States_Numbers,SPSD_Stable_State_Numbers,Num_var,Pr_Gate_Mtx)

Next_Dest_State = [];
Next_Initial_State = [];
Stable_Observed_State = [];

%Find the gate that are not still detected certainly (their probability is less than 1)
Uncertain_Index = find(Pr_Input_State_Mtx(:,1)<1);

if ~isempty(Uncertain_Index)
    m=0;
    for Current_State_Number = SPSD_Observed_States_Numbers'
        if SPSD_Stable_State_Numbers(Current_State_Number+1) %The observed state should be stable to be used as the initial state
            Stable_Observed_State(end+1) = Current_State_Number;
            m=m+1;
            Current_State = de2bi(Current_State_Number,Num_var,'left-msb');
            %Starting from the current state, find the required D and O experiment for switching to other states
            Current_State = repmat(Current_State,length(Uncertain_Index),1);
            Perturbation_Mtx = Pr_Input_State_Mtx(Uncertain_Index,2:Num_var+1)-Current_State;
            
            Perturbation_Weight_Mtx = Perturbation_Mtx;
            
            %Calculate and sort the cost of each state transition
            Perturbation_Weight_Mtx(Perturbation_Mtx==-1)= W_D;
            Perturbation_Weight_Mtx(Perturbation_Mtx==1)= W_O;
            Perturbation_Cost(:,m) = sum(Perturbation_Weight_Mtx,2);
            [Perturbation_Cost(:,m),Perturbation_Cost_indx(:,m)] = sort(Perturbation_Cost(:,m));
        end
        
    end
    
    Gate_str_sorted = {};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  More than one states are candidate as the initial state
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Perturbation_Cost,2)>1
        %Choose the best candidtae for next perturbation
        if size(Perturbation_Cost,1)>1
            Min_Perturbation_Cost = min(Perturbation_Cost);
        else
            Min_Perturbation_Cost = Perturbation_Cost;
        end
        
        [Min_Perturbation_Cost , min_indx] = sort(Min_Perturbation_Cost);
        
        %We keep every path with the minimum cost equal to the first element in the sorted costs
        min_indx(Min_Perturbation_Cost ~= Min_Perturbation_Cost(1)) = [];

        Perturbation_Cost = Perturbation_Cost(:,min_indx);
        Perturbation_Cost_indx = Perturbation_Cost_indx(:,min_indx);
        Stable_Observed_State = Stable_Observed_State(min_indx);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  More than one states are candidate as the lowest-cost destination state
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(min_indx)>1
            %Select the state with the lowest cost summation to all other destination states
            if size(Perturbation_Cost,1)>1
                Sum_Cost = sum(Perturbation_Cost);
            else
                Sum_Cost = Perturbation_Cost;
            end
            
            %TBD: If the mins are equal we should select from them.
            [~, min_sum_cost_indx] = min(Sum_Cost);
            Perturbation_Cost = Perturbation_Cost(:,min_sum_cost_indx);
            
            Next_Initial_State = de2bi(Stable_Observed_State(min_sum_cost_indx),Num_var,'left-msb');
            %Make the output table
            for i=1:size(Perturbation_Cost_indx,1)
                Index_in_Pr_Input_State_Mtx =  Uncertain_Index(Perturbation_Cost_indx(i,min_sum_cost_indx));
                Gate_sorted_cost{i,1} = cell2mat(Pr_Gate_Mtx{Index_in_Pr_Input_State_Mtx});
                if i==1
                    Next_Dest_State = Pr_Input_State_Mtx(Index_in_Pr_Input_State_Mtx,2:end);
                end
            end
            Pr_Input_State_Mtx(Uncertain_Index,:) = Pr_Input_State_Mtx(Uncertain_Index(Perturbation_Cost_indx(:,min_sum_cost_indx)),:);
            
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Only one state is candidate as the lowest-cost destination state
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Next_Initial_State = de2bi(Stable_Observed_State,Num_var,'left-msb');
            %Make the output table
            for i=1:size(Perturbation_Cost_indx,1)
                Index_in_Pr_Input_State_Mtx =  Uncertain_Index(Perturbation_Cost_indx(i));
                Gate_sorted_cost{i,1} = cell2mat(Pr_Gate_Mtx{Index_in_Pr_Input_State_Mtx});
                if i==1
                    Next_Dest_State = Pr_Input_State_Mtx(Index_in_Pr_Input_State_Mtx,2:end);
                end
            end
            Pr_Input_State_Mtx(Uncertain_Index,:) = Pr_Input_State_Mtx(Uncertain_Index(Perturbation_Cost_indx),:);
            
        end
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Only one state is candidate as the initial state
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Next_Initial_State = de2bi(Stable_Observed_State,Num_var,'left-msb');
        %Make the output table
        for i=1:size(Perturbation_Cost_indx,1)
            Index_in_Pr_Input_State_Mtx =  Uncertain_Index(Perturbation_Cost_indx(i));
            Gate_sorted_cost{i,1} = cell2mat(Pr_Gate_Mtx{Index_in_Pr_Input_State_Mtx});
            if i==1
                Next_Dest_State = Pr_Input_State_Mtx(Index_in_Pr_Input_State_Mtx,2:end);
            end
            
        end
        Pr_Input_State_Mtx(Uncertain_Index,:) = Pr_Input_State_Mtx(Uncertain_Index(Perturbation_Cost_indx),:);
        
    end
    
    Gate_Cost_Table = table(Perturbation_Cost,'RowNames',Gate_sorted_cost);
    
    
else
    
    Gate_sorted_cost{1} = cell2mat(Pr_Gate_Mtx{1});
    Perturbation_Cost = 0;
    Gate_Cost_Table = table(Perturbation_Cost,'RowNames',Gate_sorted_cost);
    
    
end

Cost = Perturbation_Cost(1);
