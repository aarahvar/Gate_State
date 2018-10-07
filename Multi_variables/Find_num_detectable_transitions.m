clc
clear States gState PIP;

%% Make the matrix of basic gates
Gate(1).gate=[1 0 0 0;0 1 1 1];
Gate(1).Name = 'AND';

Gate(2).gate=[0 1 1 1;1 0 0 0];
Gate(2).Name = 'NAND';

Gate(3).gate=[1 1 1 0;0 0 0 1];
Gate(3).Name = 'OR';

Gate(4).gate=[0 0 0 1;1 1 1 0];
Gate(4).Name = 'NOR';

Gate(5).gate=[1 0 1 1;0 1 0 0];
Gate(5).Name = 'IMPLY1';

Gate(6).gate=[0 1 0 0;1 0 1 1];
Gate(6).Name = 'NIMPLY1';

Gate(7).gate=[1 1 0 1;0 0 1 0];
Gate(7).Name = 'IMPLY2';

Gate(8).gate=[0 0 1 0;1 1 0 1];
Gate(8).Name = 'NIMPLY2';

Num_var =  5 ;

Num_Gate = 8;

Num_State = 2^Num_var;
States = zeros(Num_State,Num_Gate^(Num_var-1));
index = zeros(1,Num_Gate^(Num_var-1));

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

%% Generate the output of gate combinations


%define different orders of combination
G = [];
switch Num_var
    case 2
        for g1 = 1:Num_Gate
            G = [G g1];
        end
        
    case 3
        for g1=1:Num_Gate
            for g2=1:Num_Gate
                G = [G g1*10+g2];
            end
        end
        
    case 4
        for g1=1:Num_Gate
            for g2=1:Num_Gate
                for g3=1:Num_Gate
                    G = [G g1*100+g2*10+g3];
                end
            end
        end
        
    case 5
        for g1=1:Num_Gate
            for g2=1:Num_Gate
                for g3=1:Num_Gate
                    for g4=1:Num_Gate
                        G = [G g1*1000+g2*100+g3*10+g4];
                    end
                end
            end
        end
        
    case 6
        for g1=1:Num_Gate
            for g2=1:Num_Gate
                for g3=1:Num_Gate
                    for g4=1:Num_Gate
                        for g5=1:Num_Gate
                            G = [G g1*10000+g2*1000+g3*100+g4*10+g5];
                        end
                    end
                end
            end
        end
        
    case 7
        for g1=1:Num_Gate
            for g2=1:Num_Gate
                for g3=1:Num_Gate
                    for g4=1:Num_Gate
                        for g5=1:Num_Gate
                            for g6=1:Num_Gate
                                G = [G g1*100000+g2*10000+g3*1000+g4*100+g5*10+g6];
                            end
                        end
                    end
                end
            end
        end
        
    case 8
        for g1=1:Num_Gate
            for g2=1:Num_Gate
                for g3=1:Num_Gate
                    for g4=1:Num_Gate
                        for g5=1:Num_Gate
                            for g6=1:Num_Gate
                                for g7=1:Num_Gate
                                    G = [G g1*1000000+g2*100000+g3*10000+g4*1000+g5*100+g6*10+g7];
                                end
                            end
                        end
                    end
                end
            end
        end
        

end


%Iterative formula for calculation of outputs (G_n is the n-th gate type#)
% S_0 = eye(4)
% O_n = G_n x S_{n-1}
% S_n = kron(O_n,eye(2))

m=0;
for s=1:length(G)
    S = eye(4);
    m=m+1;
    G_str = num2str(G(s));
    for i=1:length(G_str)
        g_index(i) = str2num(G_str(i));
    end
    
    for i=1:(Num_var-1)
        O = Gate(g_index(i)).gate*S;
        S = kron(O,eye(2));
    end
    
    States(:,m)=O(1,:)';
end

switch Num_var
    case 2
         States = [States;G;zeros(3,length(index))];
    case 3
         States = [States;G;zeros(3,length(index))];
    case 4
        for s=1:length(G)
            m=m+1;
            G_str = num2str(G(s));
            for i=1:length(G_str)
                g_index(i) = str2num(G_str(i));
            end
            
            O=Gate(g_index(3)).gate*kron(Gate(g_index(1)).gate,Gate(g_index(2)).gate);
            
            States(:,m)=O(1,:)';
        end
        States = [States;[G G*10];zeros(7,2*length(index))];

    case 5
        for s=1:length(G)
            m=m+1;
            G_str = num2str(G(s));
            for i=1:length(G_str)
                g_index(i) = str2num(G_str(i));
            end
            
            O=Gate(g_index(4)).gate*kron(Gate(g_index(3)).gate*kron(Gate(g_index(1)).gate,Gate(g_index(2)).gate),eye(2));
            
            States(:,m)=O(1,:)';
        end
        
        for s=1:length(G)
            m=m+1;
            G_str = num2str(G(s));
            for i=1:length(G_str)
                g_index(i) = str2num(G_str(i));
            end
            
            O = Gate(g_index(4)).gate*kron(Gate(g_index(2)).gate*kron(Gate(g_index(1)).gate,eye(2)),Gate(g_index(3)).gate);
            
            States(:,m)=O(1,:)';
        end
        
        
        
        States = [States;[G G*10 G*100];zeros(23,3*length(index))];
        
end



%% make the gate number code of similar combination outputs
States = unique(States','rows')';


c= size(States,2);
i=1;
while i<c
    Family_Member_Index = find(sum(abs(States(1:2^Num_var,:)-repmat(States(1:2^Num_var,i),1,c)))==0)';
    States(2^Num_var+1:2^Num_var+length(Family_Member_Index),i)=States(2^Num_var+1,Family_Member_Index);
    States(:,Family_Member_Index(2:end))=[];
    i=i+1;
    c= size(States,2);
    
end



States = [States;zeros(4,size(States,2))];


%% Analyze only knockdown and knockdown+over-expression to find how many gates can be determined by only knockdown and knockdown+over-expression

% This variable show how many changes in output is visible for each state
% for only knockdowns
Num_out_change_D = zeros(Num_State,size(States,2));
Num_out_change_DO = zeros(Num_State,size(States,2));
Num_out_change_D_PPI = zeros(Num_State,size(States,2));
Num_out_change_D_O_PPI = zeros(Num_State,size(States,2));

for c=1:size(States,2)
    T = States(1:Num_State,c);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:Num_var
        D_transit(gState(i).Transit_source_D,i) = xor(T(gState(i).Transit_source_D),T(gState(i).Transit_dest_D));
    end
    
    Num_out_change_D(1:Num_State,c)=sum(D_transit,2);
    States(end-3,c)= sum(Num_out_change_D(1:Num_State,c)== Num_var);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %D+O
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:Num_var
        O_transit(gState(i).Transit_source_O,i) = xor(T(gState(i).Transit_source_O),T(gState(i).Transit_dest_O));
    end
    
    Num_out_change_DO(1:Num_State,c)=sum(O_transit,2)+sum(D_transit,2);
    States(end-2,c)= sum(Num_out_change_DO(1:Num_State,c)== Num_var);
    
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
    
    %Create input bit combinations
    Input_States = de2bi(Num_State-1:-1:0,'left-msb');
    
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
    
    States(end-1,c)= sum(Num_out_change_D_PPI(1:Num_State,c)== Num_var);
    clear D_transit_D_PPI_Temp;
    
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
            %Zer values are the inputs that the current PPI does not have any
            %effect on the revealing of them
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
    
    
    
    
    
    
    
    %     for i=1:Num_var
    %         O_transit_O_PPI(gState(i).Transit_source_O,i) = O_transit_O_PPI(gState(i).Transit_source_O,i)+ (PPI(sub2ind(size(PPI), gState(i).Transit_source_O, gState(i).Transit_dest_O)).*~O_transit_O_PPI(gState(i).Transit_source_O,i));
    %     end
    %
    %     O_transit_O_PPI = O_transit_O_PPI>0;
    %
    Num_out_change_D_O_PPI(1:Num_State,c) =sum(O_transit_O_PPI,2);
    States(end,c)= sum(Num_out_change_D_O_PPI(1:Num_State,c)== Num_var);
    clear O_transit_O_PPI_Temp;
    
end




%% Summerize results
% Row_Name = {'111';'110';'101';'100';'011';'010';'001';'000';'Comb 1';'Comb 2';'D';'D+O';'D+PPI';'D+O+PPI'};
% Table = table(States(:,1),States(:,2),States(:,3),States(:,4),States(:,5),States(:,6),States(:,7),States(:,8),States(:,9),States(:,10),States(:,11),States(:,12),...
%     States(:,13),States(:,14),States(:,15),States(:,16),States(:,17),States(:,18),States(:,19),States(:,20),States(:,21),States(:,22),States(:,23),States(:,24),States(:,25),...
%     States(:,26),States(:,27),States(:,28),States(:,29),States(:,30),States(:,31),States(:,32),'RowNames',Row_Name)
%
Row_Name = {'D';'D+O';'D+PPI';'D+O+PPI';'*+cis'};
%In each transit case, we can find only Num_Unique_State_Transition(1) states out of size(States,2) states that if
%we are there, we can determine the gate combinations uniquely
UB_Avg_P_Determination = [sum(States(end-3:end,:),2)/size(States,2)/Num_State*100;100];
UB_Avg_P_Determination = strcat(num2str(round(UB_Avg_P_Determination*100)/100),'%');
Table_Edge = table(UB_Avg_P_Determination,'RowNames',Row_Name)

