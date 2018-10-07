clc
clear States gState PIP;

load('Target_Output_Full.mat');

%%
Num_var = 3 ;

Num_Gate = 10 ;  %Num_Gate=8 ignores XOR and XNOR

Num_State = 2^Num_var;
index = zeros(1,Num_Gate^(Num_var-1));

%Create input bit combinations
Input_States = de2bi(Num_State-1:-1:0,'left-msb');


%% Define knockdown matrices
for i=0:Num_var-1
    gState(Num_var-i).D = zeros(Num_State);
    dest_indx = zeros(1,2^i);
    src_ind = zeros(1,2^i);
    while dest_indx(end)<Num_State
        src_ind =  dest_indx(end) +(1:2^i);
        dest_indx = src_ind(end)+ (1:2^i);
        
        for j=1:2^i
            gState(Num_var-i).D(src_ind(j),dest_indx(j))=1;
        end
        
    end
    
    gState(Num_var-i).Transit_source_D = find(sum(gState(Num_var-i).D,2));
    gState(Num_var-i).Transit_dest_D = find(sum(gState(Num_var-i).D))';
    
    gState(Num_var-i).O = gState(Num_var-i).D';
    gState(Num_var-i).Transit_source_O = gState(Num_var-i).Transit_dest_D;
    gState(Num_var-i).Transit_dest_O = gState(Num_var-i).Transit_source_D;
    
end


D_transit = zeros(Num_State,Num_var);

O_transit = zeros(Num_State,Num_var);

%% Define PPI matrices
%The PPI matrix shows that if we have PP intercations, it can it can
%reveale which state transition. The (imj) element in PPI matrix is 1 if a
%PPI can reveal that transition
b = dec2bin([0:2^Num_var- 1])=='0';
PPI = zeros(Num_State);
for i=1:Num_var-1
    for j=i+1:Num_var
        %Find cases that only one state is one and only one state is 1 (the
        %possible cases that PPI can provide more information)
        PPI_Index_1 = find(((b(:,i)&b(:,j))==0)&((b(:,i)|b(:,j))==1));
        
        %Here we want to find the state transitions due to PPIs. If PPI
        %reveals a transitions, it should transfer us to state 11 for those
        %inputs. Hence, we find the closes 11 state to the current state
        PPI_Full_Interact = b';
        
        %Make the state of the interacted genes to 11
        PPI_Full_Interact([i j],:)=1;
        for k=1:length(PPI_Index_1)
            %Find the closes 11 state#
            indx = find(sum(abs(PPI_Full_Interact-repmat(PPI_Full_Interact(:,PPI_Index_1(k)),1,size(PPI_Full_Interact,2))))==0);
            PPI(PPI_Index_1(k),indx(1))=1;
        end
        
    end
end

PPI = (PPI + PPI')>0;

States = Target_Output(Num_var).States;
States = [States;zeros(4,size(States,2))];




%% Find the index of the states of each combination that are unique between all combinations and by observing
Interact = 0;
for(i=1:Num_var)
    Interact = Interact+gState(i).D +  gState(i).O;
end

%those transitions we can deremine the type of the combinaation
Num_Determining_States_In_Each_Combination = zeros(Num_var,size(States,2));
%This matrix keeps the index of states of each combination family having unique transitions
Gate_Combination_Unique_Transition_indx = zeros(3,size(States,2));

%This matrix has the same size as "States" matrix and for each combination
%family (columns), the states that have unique transitions are set to 1
Gate_Combination_Unique_Transition_Matrix = zeros(Num_State,size(States,2));

ind = [];
for i=1:Num_State
    cnt = 0;
    for dbl_indx=find(Interact(i,:))
        cnt = cnt+1;
        index = [i find(Interact(i,:)) dbl_indx find(Interact(dbl_indx,:))];
        Sub_State = States(index,:);
        for j=1:size(Sub_State,2)
            %Find the state transitions that are happening only once
            num_state = sum(find(sum(abs(Sub_State-repmat(Sub_State(:,j),1,size(Sub_State,2))))==0)~=0);
            ind(cnt,j) = num_state;
            if num_state ==1
                Num_Determining_States_In_Each_Combination(cnt,j) = Num_Determining_States_In_Each_Combination(cnt,j)+1;
                Gate_Combination_Unique_Transition_indx(Num_Determining_States_In_Each_Combination(cnt,j),j)=i;
                Gate_Combination_Unique_Transition_Matrix(i,j)=1;
            end
        end
    end
end


%% Calculate the Upper Bound of probability of unique determination of different gate combination families

% These variable show how many of the edges of each state of each combination family are approved in the state diagram
Num_out_change_D = zeros(Num_State,size(States,2));
Num_out_change_DO = zeros(Num_State,size(States,2));
Num_out_change_D_PPI = zeros(Num_State,size(States,2));
Num_out_change_D_O_PPI = zeros(Num_State,size(States,2));


for c=1:size(States,2)
    T = States(1:Num_State,c);
    
    %D
    for i=1:Num_var
        %find the changes in outputs for each state of the combination family
        D_transit(gState(i).Transit_source_D,i) = xor(T(gState(i).Transit_source_D),T(gState(i).Transit_dest_D));
    end
    
    %count the number of changes in the output of each state
    Num_out_change_D(1:Num_State,c)=sum(D_transit,2);
    %Find the states that they have changes on all edges while those states are determining states
    States(end-3,c)= sum(Num_out_change_D(1:Num_State,c).*Gate_Combination_Unique_Transition_Matrix(:,c)== Num_var);
    
    
    % D+PPI
    D_transit_D_PPI = D_transit;
    for i=1:Num_var
        %We detect transitions that are due to D OR PPI
        %D_transit_D_PPI(gState(i).Transit_source_D,i) = D_transit_D_PPI(gState(i).Transit_source_D,i)+ (PPI(sub2ind(size(PPI), gState(i).Transit_source_D, gState(i).Transit_dest_D)).*~D_transit_D_PPI(gState(i).Transit_source_D,i));
        D_transit_D_PPI(gState(i).Transit_source_D,i) = D_transit_D_PPI(gState(i).Transit_source_D,i) | (PPI(sub2ind(size(PPI), gState(i).Transit_source_D, gState(i).Transit_dest_D)));
    end
    
    %count the number of changes in the output of each state
    Num_out_change_D_PPI(1:Num_State,c) =sum(D_transit_D_PPI,2);
    %Find the states that they have changes on all edges while those states are determining states
    States(end-2,c)= sum(Num_out_change_D_PPI(1:Num_State,c).*Gate_Combination_Unique_Transition_Matrix(:,c)== Num_var);
    %D+O
    for i=1:Num_var
        %find the changes in outputs for each state of the combination family
        O_transit(gState(i).Transit_source_O,i) = xor(T(gState(i).Transit_source_O),T(gState(i).Transit_dest_O));
    end
    
    %count the number of changes in the output of each state
    Num_out_change_DO(1:Num_State,c)=sum(O_transit,2)+sum(D_transit,2);
    %Find the states that they have changes on all edges while those states are determining states
    States(end-1,c)= sum(Num_out_change_DO(1:Num_State,c).*Gate_Combination_Unique_Transition_Matrix(:,c)== Num_var);
    
    % D+O+PPI
    O_transit_O_PPI = O_transit+D_transit_D_PPI;
    for i=1:Num_var
        %Here, we consider interactions in PPI matrix that their corresponding edge has not been considered in "D+O analysis"
        O_transit_O_PPI(gState(i).Transit_source_O,i) = O_transit_O_PPI(gState(i).Transit_source_O,i)+ (PPI(sub2ind(size(PPI), gState(i).Transit_source_O, gState(i).Transit_dest_O)).*~O_transit_O_PPI(gState(i).Transit_source_O,i));
    end
    
    %The elements on O_transit_O_PPI matrix should be only 1 or 0
    O_transit_O_PPI = O_transit_O_PPI>0;
    
    %count the number of changes in the output of each state
    Num_out_change_D_O_PPI(1:Num_State,c) =sum(O_transit_O_PPI,2);
    %Find the states that they have changes on all edges while those states are determining states
    States(end,c)= sum(Num_out_change_D_O_PPI(1:Num_State,c).*Gate_Combination_Unique_Transition_Matrix(:,c)== Num_var);
    
end




%% Summerize results
% Row_Name = {'111';'110';'101';'100';'011';'010';'001';'000';'Comb 1';'Comb 2';'D';'D+O';'D+PPI';'D+O+PPI'};
% Table = table(States(:,1),States(:,2),States(:,3),States(:,4),States(:,5),States(:,6),States(:,7),States(:,8),States(:,9),States(:,10),States(:,11),States(:,12),...
%     States(:,13),States(:,14),States(:,15),States(:,16),States(:,17),States(:,18),States(:,19),States(:,20),States(:,21),States(:,22),States(:,23),States(:,24),States(:,25),...
%     States(:,26),States(:,27),States(:,28),States(:,29),States(:,30),States(:,31),States(:,32),'RowNames',Row_Name)
%
Row_Name = {'D';'D+PPI';'D+O';'D+O+PPI';'D+O+cis'};
%In each transit case, we can find only Num_Unique_State_Transition(1) states out of size(States,2) states that if
%we are there, we can determine the gate combinations uniquely
UB_Avg_P_Determination = [sum(States(end-3:end,:),2)/size(States,2)/Num_State*100;100*mean(Num_Determining_States_In_Each_Combination/Num_State)];
UB_Avg_P_Determination = strcat(num2str(round(UB_Avg_P_Determination*100)/100),'%');
Table_Edge = table(UB_Avg_P_Determination,'RowNames',Row_Name)

