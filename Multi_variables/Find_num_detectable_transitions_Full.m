clc
clear States gState PIP ;

load('Target_Output_Full.mat');
%% Make the matrix of basic gates


Num_var = 2 ;

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

load('Target_Output_Full.mat');
States = Target_Output(Num_var).States;

States = [States;zeros(4,size(States,2))];



%% Analyze only knockdown and knockdown+over-expression to find how many gates can be determined by only knockdown and knockdown+over-expression

% This variable show how many changes in output is visible for each state
Num_out_change_D = zeros(Num_State,size(States,2));
Num_out_change_DO = zeros(Num_State,size(States,2));
Num_out_change_D_PPI = zeros(Num_State,size(States,2));
Num_out_change_D_O_PPI = zeros(Num_State,size(States,2));

%The following structure shows the state number of each target output that
%we have evidence for all TFs for different scenarios. For example:
%
% Visible_TF_Effect_Index(3).D(1).Num_Output: shows the order of output
% valuse (column #) in Target_Output_Full
%
%Visible_TF_Effect_Index(3).D(1).Index: Shows the indexes of the above
%output that they show evidences for all TFs (i.e., if we are in other
%states not belonging to this set, we cannot provide enough evidence for
%existence of all TFs)
Visible_TF_Effect_Index(Num_var) = struct('D',[],'DO',[],'D_PPI',[],'DO_PPI',[],'cis',[]);

    for i=1:Num_State
        Visible_TF_Effect_Index(Num_var).D(i).Output_index=[];
        Visible_TF_Effect_Index(Num_var).DO(i).Output_index=[];
        Visible_TF_Effect_Index(Num_var).D_PPI(i).Output_index=[];
        Visible_TF_Effect_Index(Num_var).DO_PPI(i).Output_index=[];
    end

for c=1:size(States,2)
    T = States(1:Num_State,c);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:Num_var
        D_transit(gState(i).Transit_source_D,i) = xor(T(gState(i).Transit_source_D),T(gState(i).Transit_dest_D));
    end
    
    
    Num_out_change_D(1:Num_State,c)=sum(D_transit,2);
    States(end-3,c)= sum(Num_out_change_D(:,c)== Num_var);
    
    if any(States(end-3,c))
        full_transition_index = find(Num_out_change_D(:,c)== Num_var)';
        for indx=full_transition_index
            Visible_TF_Effect_Index(Num_var).D(indx).Output_index(end+1) = c;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %D+O
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:Num_var
        O_transit(gState(i).Transit_source_O,i) = xor(T(gState(i).Transit_source_O),T(gState(i).Transit_dest_O));
    end
    
    Num_out_change_DO(1:Num_State,c)=sum(O_transit,2)+sum(D_transit,2);
    States(end-1,c)= sum(Num_out_change_DO(1:Num_State,c)== Num_var);

    if any(States(end-1,c))
        full_transition_index = find(Num_out_change_DO(:,c)== Num_var)';
        for indx=full_transition_index
            Visible_TF_Effect_Index(Num_var).DO(indx).Output_index(end+1) = c;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % D+PPI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %To measure the effect of PPI on the probability of deection of TFs, we
    %consider the PPI as a revealer PPI if the protein of this PPI attaches
    %to one of the previousely detected D transitions that make an output
    %change. If we have some PPIs without any evidence about the attachment
    %of them to a known TF, we do not consider it as an effective TF.
    
    
    %Since PPI matix is symmetric, we just consider the upper side
    PPI_Ind = find(triu(PPI));
    
    %Find the row(R) and column(C) index of 1s in PPI matrix
    [R,C]=ind2sub(size(PPI),PPI_Ind);
    
    
    %D_transit is a Num_State x Num_var matrix that shows the status of the
    %detection of each Di for that state saying that if it is making an
    %output change or not
    %In D_transit_D_PPI we will have both detected transitions by D and PPI
    D_transit_D_PPI = D_transit;
    D_transit_D_PPI_Temp = zeros(size(D_transit_D_PPI));
    
    %The following While loop is for the cases that after the first
    %iteration and revealing some TFs, the new ones may reveal some others
    %in the next iteration. We continue until we do not see any change in
    %the result (Running the code showed that more iteration (for D+PPI) does not add
    %any new TF or transition)
    m=0;
    while(sum(sum(D_transit_D_PPI_Temp-D_transit_D_PPI))~=0)
        m=m+1;
        
        D_transit_D_PPI_Temp = D_transit_D_PPI;
        
        
        for i=1:length(PPI_Ind)
            
            
            %Find the input bit that can be revealed by PPI (it should be
            %toggled between 0 and 1)
            Candidate_PPI = xor(Input_States(R(i),:) , Input_States(C(i),:));
            Candidate_PPI_Ind = find(Candidate_PPI);
            %Zer values are the inputs that the current PPI does not have any
            %effect on the revealing of them
            Non_Candidate_PPI_Ind = find(Candidate_PPI==0);
            
            if D_transit_D_PPI(R(i),Candidate_PPI_Ind)==1
                continue;
            end
            
            %If for any potential candidate of revealing by PPI the the other
            %partner of PPI has been detected before, we will consider this one
            %as a detectable transition by PPI
            
            for j=1:length(Non_Candidate_PPI_Ind)
                if(D_transit_D_PPI(R(i),Non_Candidate_PPI_Ind(j))==1)
                    D_transit_D_PPI(R(i),Candidate_PPI_Ind)=1;
                    break;
                end
            end
        end
    end
    
    Num_out_change_D_PPI(1:Num_State,c) =sum(D_transit_D_PPI,2);
    
    States(end-2,c)= sum(Num_out_change_D_PPI(1:Num_State,c)== Num_var);
 
    if any(States(end-2,c))
        full_transition_index = find(Num_out_change_D_PPI(:,c)== Num_var)';
        for indx=full_transition_index
            Visible_TF_Effect_Index(Num_var).D_PPI(indx).Output_index(end+1) = c;
        end
    end
    
    
    clear D_transit_D_PPI_Temp  ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % D+O+PPI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    O_transit_O_PPI = O_transit+D_transit;
    O_transit_O_PPI_Temp = zeros(size(O_transit_O_PPI));
    
    m=0;
    while(sum(sum(O_transit_O_PPI_Temp-O_transit_O_PPI))~=0)
        m=m+1;
        
        O_transit_O_PPI_Temp = O_transit_O_PPI;
        
        
        for i=1:length(PPI_Ind)
            
            
            %Find the input bit that can be revealed by PPI (it should be
            %toggled between 0 and 1)
            Candidate_PPI = xor(Input_States(R(i),:) , Input_States(C(i),:));
            Candidate_PPI_Ind = find(Candidate_PPI);
            %Zero values are the inputs that the current PPI does not have any
            %effect on revealing of them
            Non_Candidate_PPI_Ind = find(Candidate_PPI==0);
            
            if O_transit_O_PPI(R(i),Candidate_PPI_Ind)==1
                continue;
            end
            
            %If for any potential candidate of revealing by PPI the the other
            %partner of PPI has been detected before, we will consider this one
            %as a detectable transition by PPI
            
            for j=1:length(Non_Candidate_PPI_Ind)
                if(O_transit_O_PPI(R(i),Non_Candidate_PPI_Ind(j))==1)
                    O_transit_O_PPI(R(i),Candidate_PPI_Ind)=1;
                    break;
                end
            end
        end
    end
    
 
    Num_out_change_D_O_PPI(1:Num_State,c) =sum(O_transit_O_PPI,2);
    States(end,c)= sum(Num_out_change_D_O_PPI(1:Num_State,c)== Num_var);

    if any(States(end,c))
        full_transition_index = find(Num_out_change_D_O_PPI(:,c)== Num_var)';
        for indx=full_transition_index
            Visible_TF_Effect_Index(Num_var).DO_PPI(indx).Output_index(end+1) = c;
        end
    end
    
    clear O_transit_O_PPI_Temp   ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          *+cis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Num_State
    Visible_TF_Effect_Index(Num_var).cis(i).Output_index = 1:size(States,2);
end


%% Summerize results
% Row_Name = {'111';'110';'101';'100';'011';'010';'001';'000';'Comb 1';'Comb 2';'D';'D+O';'D+PPI';'D+O+PPI'};
% Table = table(States(:,1),States(:,2),States(:,3),States(:,4),States(:,5),States(:,6),States(:,7),States(:,8),States(:,9),States(:,10),States(:,11),States(:,12),...
%     States(:,13),States(:,14),States(:,15),States(:,16),States(:,17),States(:,18),States(:,19),States(:,20),States(:,21),States(:,22),States(:,23),States(:,24),States(:,25),...
%     States(:,26),States(:,27),States(:,28),States(:,29),States(:,30),States(:,31),States(:,32),'RowNames',Row_Name)
%
Row_Name = {'D';'D+PPI';'D+O';'D+O+PPI';'*+cis'};
%In each transit case, we can find only Num_Unique_State_Transition(1) states out of size(States,2) states that if
%we are there, we can determine the gate combinations uniquely
UB_Avg_P_Determination = [sum(States(end-3:end,:),2)/size(States,2)/Num_State*100;100];
UB_Avg_P_Determination = strcat(num2str(round(UB_Avg_P_Determination*100)/100),'%');
Table_Edge = table(UB_Avg_P_Determination,'RowNames',Row_Name)

