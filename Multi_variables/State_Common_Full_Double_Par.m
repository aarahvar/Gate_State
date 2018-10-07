clc
clear States gState PIP Num_common_States index Mean_Common_State;

load('Target_Output_Full.mat');
tic
%%
Num_var = 2 ;

Num_Gate = 10 ;  %Num_Gate=8 ignores XOR and XNOR

Num_State = 2^Num_var;

%Create input bit combinations
Input_States = de2bi(Num_State-1:-1:0,'left-msb');

Data_Type_Array = {'D','DO'};
Data_Type = Data_Type_Array{1};


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
switch Data_Type
    case 'D'
        for(i=1:Num_var)
            Interact = Interact+gState(i).D;
        end
        
    case 'DO'
        for(i=1:Num_var)
            Interact = Interact+gState(i).D +  gState(i).O;
        end
end

% % %those transitions we can deremine the type of the combinaation
% % Num_Determining_States_In_Each_Combination = zeros(1,size(States,2));
% % %This matrix keeps the index of states of each combination family having unique transitions
% % Gate_Combination_Unique_Transition_indx = zeros(3,size(States,2));
% % 
% % %This matrix has the same size as "States" matrix and for each combination
% % %family (columns), the states that have unique transitions are set to 1
% % Gate_Combination_Unique_Transition_Matrix = zeros(Num_State,size(States,2));

for i=1:Num_State
    Num_common_States(i).Num_State = zeros(1,size(States,2));
end

cnt = 0;
for i=1:Num_State
    if isempty(find(Interact(i,:)))
        cnt = cnt+1;
        index(cnt).index = i;
        continue;
    end
    for dbl_indx=find(Interact(i,:))
        cnt = cnt+1;
        if isempty(find(Interact(dbl_indx,:)))
            index(cnt).index = [i find(Interact(i,:))];
            continue;
        end
        index(cnt).index = [i find(Interact(i,:)) dbl_indx find(Interact(dbl_indx,:)) ];
    end
    Mean_Common_State(i) = mean(min(Num_common_States(i).Num_State));
end


num_state = zeros(1,size(States,2));
for i=1:length(index)
    Sub_State = States(unique(index(i).index),:);
    parfor j=1:size(Sub_State,2)
        %Find the state transitions that are happening only once
        if size(Sub_State,1)==1
            num_state(j) = sum(find((abs(Sub_State-repmat(Sub_State(:,j),1,size(Sub_State,2))))==0)~=0);
        else
            num_state(j) = sum(find(sum(abs(Sub_State-repmat(Sub_State(:,j),1,size(Sub_State,2))))==0)~=0);
        end
        
    end
    Num_common_States(index(i).index(1)).Num_State(i,:) = num_state;
end

for i=1:Num_State
    A = Num_common_States(i).Num_State;
    A(A==0)=NaN;
    if size(A,1)==1
        Min_A = A;
    else
        Min_A = min(A);
    end
    
    Min_A(isnan(Min_A))=0;
    Mean_Common_State(i) = mean(Min_A);

end

%% Calculate the Upper Bound of probability of unique determination of different gate combination families
Mean_Common_State
mean(Mean_Common_State)
toc
