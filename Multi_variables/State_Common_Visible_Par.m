

clc
clear all
load('Target_Output_Full.mat');
load('Visible_TF_Effect_Index.mat');
%% Make the matrix of basic gates

Num_var =  3 ;

Num_State = 2^Num_var;

Input_States = de2bi(Num_State-1:-1:0,'left-msb');

Data_Type_Array = {'D','D_PPI','D_cis','DO','DO_PPI','DO_cis'};
Data_Type = Data_Type_Array{4};

%% Generate the output of gate combinations

States = Target_Output(Num_var).States;
clear Target_Output;

States = [States;zeros(4,size(States,2))];

%% Define knockdown matrices
%If in the state diagram, there is a directed edge from node i to node j,
%we put i in Transit_source_D or Transit_source_O and put j in
%Transit_dest_D or Transit_dest_O
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
        %inputs. Hence, we find the closest 11 state to the current state
        PPI_Full_Interact = b';
        
        %Make the state of the interacted genes to 11
        PPI_Full_Interact([i j],:)=1;
        for k=1:length(PPI_Index_1)
            %Find the closes 11 state#
            idx = find(sum(abs(PPI_Full_Interact-repmat(PPI_Full_Interact(:,PPI_Index_1(k)),1,size(PPI_Full_Interact,2))))==0);
            PPI(PPI_Index_1(k),idx(1))=1;
        end
        
    end
end

PPI = (PPI + PPI')>0;




%% Find the index of the states of each combination that are unique between all combinations and by observing
Interact = 0;

switch Data_Type
    case {'D','D_PPI','D_cis'}
        for(i=1:Num_var)
            Interact = Interact+gState(i).D;
        end
        
    otherwise
        for(i=1:Num_var)
            Interact = Interact+gState(i).D +  gState(i).O;
        end
end





for i=1:Num_State
    Num_common_States(i).Num_State = [];
end

cnt = 0;
for i=1:Num_State
    index(i).index = [i find(Interact(i,:))];;
end


num_state = zeros(1,size(States,2));
for i=1:length(index)
    Sub_State = Creat_Sub_State(States,Visible_TF_Effect_Index(Num_var),Data_Type,unique(index(i).index),index(i).index(1));
    for j=1:size(Sub_State,2)
        %Find the state transitions that are happening only once
        if size(Sub_State,1)==1
            num_state(j) = sum(find((abs(Sub_State-repmat(Sub_State(:,j),1,size(Sub_State,2))))==0)~=0);
        else
            num_state(j) = sum(find(sum(abs(Sub_State-repmat(Sub_State(:,j),1,size(Sub_State,2))))==0)~=0);
        end
    end
    Num_common_States(index(i).index(1)).Num_State(end+1,:) = num_state(1:size(Sub_State,2));
end

for i=1:Num_State
    Num_Common = Num_common_States(i).Num_State;
    if isempty(Num_Common)
        Mean_Common_State(i) = 0;
        Sum_P_Detection_in_Visibles(i)=0;
        continue;
    end
    if size(Num_Common,1)==1
        Min_Num_Common = Num_Common;
    else
        Min_Num_Common = min(Num_Common);
    end
    
    %Average # of Boolean functions having common transitions for each state
    Mean_Common_State(i) = mean(Min_Num_Common);
    
    %
    Sum_P_Detection_in_Visibles(i) = sum(1./Min_Num_Common);

end

% Number of Boolean functions whose s-th state is fully visible
Num_State_Members = Get_Num_State_Members(Visible_TF_Effect_Index(Num_var),Data_Type,Num_State);


% Calculate the Upper Bound of probability of unique determination of different gate combination families
Mean_Common_State
Sum_P_Detection_in_Visibles
% Num_State_Members
P_Visibility_in_State = Num_State_Members/size(States,2)
ind = (Sum_P_Detection_in_Visibles~=0);
Unambiguity_Score = sum(P_Visibility_in_State(ind).*Sum_P_Detection_in_Visibles(ind)./Num_State_Members(ind))/Num_State
Mean_Common_in_Visible = sum(P_Visibility_in_State(ind)/sum(P_Visibility_in_State(ind)).*Mean_Common_State(ind))

return


%%
load('Common_State_Visible.mat')
Common_State_Visible(Num_var).Single_DO_cis.Unambiguity_Score = Unambiguity_Score;
Common_State_Visible(Num_var).Single_DO_cis.Mean_Common_State = Mean_Common_State;
Common_State_Visible(Num_var).Single_DO_cis.Num_State_Members = Num_State_Members;
Common_State_Visible(Num_var).Single_DO_cis.index = index;
Common_State_Visible(Num_var).Single_DO_cis.Mean_Common_in_Visible = Mean_Common_in_Visible;
Common_State_Visible(Num_var).Single_DO_cis.Sum_P_Detection_in_Visibles = Sum_P_Detection_in_Visibles;
if isfield(Common_State_Visible(Num_var).Single_DO_cis,'Visibility_Unambiguity_Score')
    Common_State_Visible(Num_var).Single_DO_cis= rmfield(Common_State_Visible(Num_var).Single_DO_cis,'Visibility_Unambiguity_Score');
end
save('Common_State_Visible.mat','Common_State_Visible')