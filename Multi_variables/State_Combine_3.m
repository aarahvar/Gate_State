clc
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

Num_var = 3;
Num_State = 2^Num_var;
States = zeros(Num_State,8^(Num_var-1));
index = zeros(1,8^(Num_var-1));

%% Define knockdown matrices
D1 = zeros(Num_State);
D2 = zeros(Num_State);
D3 = zeros(Num_State);

D1(1,5)=1;
D1(3,7)=1;
D1(2,6)=1;
D1(4,8)=1;

D2(1,3)=1;
D2(5,7)=1;
D2(2,4)=1;
D2(6,8)=1;

D3(1,2)=1;
D3(3,4)=1;
D3(5,6)=1;
D3(7,8)=1;

D1_Transit_source = find(sum(D1,2));
D1_Transit_dest = find(sum(D1))';

D2_Transit_source = find(sum(D2,2));
D2_Transit_dest = find(sum(D2))';

D3_Transit_source = find(sum(D3,2));
D3_Transit_dest = find(sum(D3))';

D_transit = zeros(Num_State,3);

%% Define over-expression matrices
O1 = D1';
O2 = D2';
O3 = D3';

O1_Transit_source = find(sum(O1,2));
O1_Transit_dest = find(sum(O1))';

O2_Transit_source = find(sum(O2,2));
O2_Transit_dest = find(sum(O2))';

O3_Transit_source = find(sum(O3,2));
O3_Transit_dest = find(sum(O3))';

O_transit = zeros(Num_State,3);

%% Define PPI matrices
PPI12 = zeros(Num_State);
PPI13 = zeros(Num_State);
PPI23 = zeros(Num_State);

PPI12(2,6)=1;
PPI12(2,4)=1;
PPI12(1,5)=1;
PPI12(1,3)=1;
PPI12 = PPI12 + PPI12';

PPI13(3,4)=1;
PPI13(3,7)=1;
PPI13(1,2)=1;
PPI13(1,5)=1;
PPI13 = PPI13 + PPI13';

PPI23(5,6)=1;
PPI23(5,7)=1;
PPI23(1,2)=1;
PPI23(1,3)=1;
PPI23 = PPI23 + PPI23';

PPI = (PPI12+PPI13+PPI23)>0;




%% Generate the output of gate combinations
for i=1:8
   Gate1_x3(i).gate= kron(Gate(i).gate,eye(2));
end

m=0;
for g1 = 1:8
    for g2 = 1:8
        m =m+1;
        A=Gate(g2).gate*Gate1_x3(g1).gate;
        States(:,m)=A(1,:)';
        index(m) = g1*10+g2;
    end
end

%% make the gate number code of similar combination outputs
States = [States;index;index];
States = unique(States','rows')';
States(end,1:2:end)=States(end,2:2:end);
States(:,2:2:end)=[];

States = [States;zeros(4,size(States,2))];



%% Analyze only knockdown and knockdown+over-expression to find how many gates can be determined by only knockdown and knockdown+over-expression

% This variable show for each state how many changes in output is visible
% for only knockdowns
Num_out_change_D = zeros(Num_State,size(States,2));
Num_out_change_DO = zeros(Num_State,size(States,2));
Num_out_change_D_PPI = zeros(Num_State,size(States,2));
Num_out_change_D_O_PPI = zeros(Num_State,size(States,2));

for c=1:size(States,2)
    T = States(1:Num_State,c);
    
    %D
    D_transit(D1_Transit_source,1) = mod(T(D1_Transit_source)+T(D1_Transit_dest),2);
    D_transit(D2_Transit_source,2) = mod(T(D2_Transit_source)+T(D2_Transit_dest),2);
    D_transit(D3_Transit_source,3) = mod(T(D3_Transit_source)+T(D3_Transit_dest),2);

    Num_out_change_D(1:Num_State,c)=sum(D_transit,2);
    States(end-3,c)= sum(Num_out_change_D(1:Num_State,c)== Num_var);
    
    %D+O
    O_transit(O1_Transit_source,1) = mod(T(O1_Transit_source)+T(O1_Transit_dest),2);
    O_transit(O2_Transit_source,2) = mod(T(O2_Transit_source)+T(O2_Transit_dest),2);
    O_transit(O3_Transit_source,3) = mod(T(O3_Transit_source)+T(O3_Transit_dest),2);
    
    Num_out_change_DO(1:Num_State,c)=sum(O_transit,2)+sum(D_transit,2);
    States(end-2,c)= sum(Num_out_change_DO(1:Num_State,c)== Num_var);
    
    % D+PPI
    D_transit_D_PPI = D_transit;
    D_transit_D_PPI(D1_Transit_source,1) = D_transit_D_PPI(D1_Transit_source,1)+ (PPI(sub2ind(size(PPI), D1_Transit_source, D1_Transit_dest)).*~D_transit_D_PPI(D1_Transit_source,1));
    D_transit_D_PPI(D2_Transit_source,2) = D_transit_D_PPI(D2_Transit_source,2)+ (PPI(sub2ind(size(PPI), D2_Transit_source, D2_Transit_dest)).*~D_transit_D_PPI(D2_Transit_source,2));
    D_transit_D_PPI(D3_Transit_source,3) = D_transit_D_PPI(D3_Transit_source,3)+ (PPI(sub2ind(size(PPI), D3_Transit_source, D3_Transit_dest)).*~D_transit_D_PPI(D3_Transit_source,3));
    
    Num_out_change_D_PPI(1:Num_State,c) =sum(D_transit_D_PPI,2);
    States(end-1,c)= sum(Num_out_change_D_PPI(1:Num_State,c)== Num_var);
    
    % D+O+PPI
    O_transit_O_PPI = O_transit+D_transit_D_PPI;
    O_transit_O_PPI(O1_Transit_source,1) = O_transit_O_PPI(O1_Transit_source,1)+ (PPI(sub2ind(size(PPI), O1_Transit_source, O1_Transit_dest)).*~O_transit_O_PPI(O1_Transit_source,1));
    O_transit_O_PPI(O2_Transit_source,2) = O_transit_O_PPI(O2_Transit_source,2)+ (PPI(sub2ind(size(PPI), O2_Transit_source, O2_Transit_dest)).*~O_transit_O_PPI(O2_Transit_source,2));
    O_transit_O_PPI(O3_Transit_source,3) = O_transit_O_PPI(O3_Transit_source,3)+ (PPI(sub2ind(size(PPI), O3_Transit_source, O3_Transit_dest)).*~O_transit_O_PPI(O3_Transit_source,3));
    O_transit_O_PPI = O_transit_O_PPI>0;
    
    Num_out_change_D_O_PPI(1:Num_State,c) =sum(O_transit_O_PPI,2);
    States(end,c)= sum(Num_out_change_D_O_PPI(1:Num_State,c)== Num_var);

end

%% Find unique states determining the gate combination type
Interact = D1+D2+D3+O1+O2+O3;
Num_Unique_State_Transition = zeros(1,Num_State);

for i=1:Num_State
    index = [i find(Interact(i,:))];
    Sub_State = States(index,:);
    Unique_Sub_State = unique(Sub_State','rows')';
    for j=1:size(Unique_Sub_State,2)
        %Find the state transitions that are happening only once
        num_state = sum(find(sum(abs(Sub_State-repmat(Unique_Sub_State(:,j),1,size(Sub_State,2))))==0)~=0);
        if num_state ==1
            Num_Unique_State_Transition(i) = Num_Unique_State_Transition(i)+1;
        end
    end
end
Num_Unique_State_Transition


%% Summerize results
Row_Name = {'111';'110';'101';'100';'011';'010';'001';'000';'Comb 1';'Comb 2';'D';'D+O';'D+PPI';'D+O+PPI'};
Table = table(States(:,1),States(:,2),States(:,3),States(:,4),States(:,5),States(:,6),States(:,7),States(:,8),States(:,9),States(:,10),States(:,11),States(:,12),...
    States(:,13),States(:,14),States(:,15),States(:,16),States(:,17),States(:,18),States(:,19),States(:,20),States(:,21),States(:,22),States(:,23),States(:,24),States(:,25),...
    States(:,26),States(:,27),States(:,28),States(:,29),States(:,30),States(:,31),States(:,32),'RowNames',Row_Name)

Row_Name = {'D';'D+O';'D+PPI';'D+O+PPI';'*+cis'};
%In each transit case, we can find only 12 states out of 32 states that if
%we are there, we can determine the gate combinations uniquely
UB_Avg_P_Determination = [sum(States(end-3:end,:),2)/32/8*(12/32)*100;100*12/32];
UB_Avg_P_Determination = strcat(num2str(round(UB_Avg_P_Determination*10)/10),'%');
Table_Edge = table(UB_Avg_P_Determination,'RowNames',Row_Name)

