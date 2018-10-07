
clc
clear SPSD;
load('Removed_Function_Index.mat');
%% Find the bit combinations (states) of input bits
%input simulation

%% ABC
SPSD(1).Initial_State_str = '011';
SPSD(1).Target_State_str = '0';

SPSD(1).Perturbed_Input_index_str = '';



%Give the Boolean function index to get the output states (Target_Output[1] and Target_Output[end] refere to 00...0 and 11..1, respectively.
Boolean_Function_Index = 22;

%This variable shows which input state combinations are stable (could be used as the initial state for transitions)
SPSD_Stable_State_Numbers = [1 1 1 1 1 1 1 1 1];


%% ABC
% SPSD(1).Initial_State_str = '010';
% SPSD(1).Target_State_str = '0';
% SPSD(1).Perturbed_Input_index_str = '1 2 3';
% SPSD(1).Output_States_After_Pert_str = '000';
%
% SPSD(2).Initial_State_str = '000';
% SPSD(2).Perturbed_Input_index_str = '1 3';
% SPSD(2).Output_States_After_Pert_str = '00';
%
% SPSD(3).Initial_State_str = '100';
% SPSD(3).Perturbed_Input_index_str = '3';
% SPSD(3).Output_States_After_Pert_str = '0';
%
% SPSD(4).Initial_State_str = '101';
% SPSD(4).Perturbed_Input_index_str = '2';
% SPSD(4).Output_States_After_Pert_str = '1';


% %% F22
% SPSD(1).Initial_State_str = '011';
% SPSD(1).Target_State_str = '1';
% SPSD(1).Perturbed_Input_index_str = '1 2 3';
% SPSD(1).Output_States_After_Pert_str = '000';
%
% SPSD(2).Initial_State_str = '111';
% SPSD(2).Perturbed_Input_index_str = '2 3';
% SPSD(2).Output_States_After_Pert_str = '11';
%
% SPSD(3).Initial_State_str = '101';
% SPSD(3).Perturbed_Input_index_str = '3';
% SPSD(3).Output_States_After_Pert_str = '0';
%
% SPSD(4).Initial_State_str = '100';
% SPSD(4).Perturbed_Input_index_str = '1';
% SPSD(4).Output_States_After_Pert_str = '0';



%% ABCD

% SPSD(1).Initial_State_str = '1111';
% %SPSD(1).Target_State_str = '1';
% SPSD(1).Perturbed_Input_index_str = '';
% %SPSD(1).Output_States_After_Pert_str = '0000';
% 
% 
% %Give the Boolean function index to get the output states (Target_Output[1] and Target_Output[end] refere to 00...0 and 11..1, respectively.
% Boolean_Function_Index = 5;
% 
% 
%  SPSD_Stable_State_Numbers = [0 0 0 0 1 1 1 1 0 0 0 0 0 1 1 1 1];
% %SPSD_Stable_State_Numbers = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
% 

%% Initialization
Num_var = length(SPSD(1).Initial_State_str) ;
Num_State = 2^Num_var;

Output_States = de2bi(Boolean_Function_Index,Num_State,'left-msb');

Gate(2).gate = {'~A~B'; '~AB';'A~B';'AB'};
Gate(3).gate = {'~A~B~C';'~A~BC';'~AB~C';'~ABC';'A~B~C';'A~BC';'AB~C';'ABC'};
Gate(4).gate = {'~A~B~C~D';'~A~B~CD';'~A~BC~D';'~A~BCD';'~AB~C~D';'~AB~CD';'~ABC~D';'~ABCD';'A~B~C~D';'A~B~CD';'A~BC~D';'A~BCD';'AB~C~D';'AB~CD';'ABC~D';'ABCD'};



%% Get Data from user

SPSD(1).Initial_State_Binary = SPSD(1).Initial_State_str=='1';

SPSD(1).Perturbed_Input_index = str2num(SPSD(1).Perturbed_Input_index_str);

SPSD(1).SPSD_Observed_States = repmat(SPSD(1).Initial_State_Binary,length(SPSD(1).Perturbed_Input_index)+1,1);
for i = 1:length(SPSD(1).Perturbed_Input_index)
    indx = SPSD(1).Perturbed_Input_index(i);
    SPSD(1).SPSD_Observed_States(i+1,indx) = ~SPSD(1).SPSD_Observed_States(i+1,indx);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                First Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate the truth table of Boolean function having common transitions %(1st step)
SPSD_Observed_States_Numbers = bi2de(SPSD(1).SPSD_Observed_States,'left-msb');
SPSD(1).Output_States = Output_States(SPSD_Observed_States_Numbers+1)';

%[SPSD(1).Boolean_Function_Common , SPSD(1).SPSD_Observed_States_Numbers] = Generate_Boolean_Function_Common_1st_Step(SPSD_Observed_States_Numbers,SPSD(1).Output_States,Num_State,Removed_Function_Index(Num_var).index);
[SPSD(1).Boolean_Function_Common , SPSD(1).SPSD_Observed_States_Numbers] = Generate_Boolean_Function_Common_1st_Step(SPSD_Observed_States_Numbers,SPSD(1).Output_States,Num_State,[]);

%Find the probability of each gate by counting the number of the ones in truth table for each input state
[Pr_Input_State_Mtx , Pr_Gate_Mtx] = Find_Gate_Probability(SPSD(1).Boolean_Function_Common,Gate,Num_var,1);

%Set the cost weight of knock-down (D) and over-expression (O)
W_D = 1;
W_O = 2;

%Find the cost of each transition from observed states as the initial state (the states who are stable)
[Gate_Cost_Table,Pr_Input_State_Mtx,Next_Initial_State,Next_Dest_State,Cost] = Find_Minimum_Cost_Transition(Pr_Input_State_Mtx,W_D,W_O,SPSD(1).SPSD_Observed_States_Numbers,SPSD_Stable_State_Numbers,Num_var,Pr_Gate_Mtx);
Total_Cost = Cost;


display(Gate_Cost_Table)
    display('-------------------------------')
    fprintf(['<strong> ' num2str(Next_Initial_State) '   ->   ' num2str(Next_Dest_State) '</strong>\n'])
    display('-------------------------------')

%%

for step = 2:Num_State
    
    SPSD(step).Initial_State_Binary = Next_Dest_State;
    %     SPSD(step).Perturbed_Input_index = str2num(SPSD(step).Perturbed_Input_index_str);
    
    
    %Generate the truth table of Boolean function having common transitions %(2nd step)
    SPSD = Generate_Boolean_Function_Common_2(SPSD,step,Output_States);
    
    %Find the probability of each gate by counting the number of the ones in truth table for each input state
    [Pr_Input_State_Mtx , Pr_Gate_Mtx] = Find_Gate_Probability(SPSD(step).Boolean_Function_Common,Gate,Num_var,step);
    
    if sum(Pr_Input_State_Mtx(:,1))==size(Pr_Input_State_Mtx,1)
        display('---------------------------')
        display('       FINISHED')
        display('---------------------------')
        display(['<strong>Total cost: ' num2str(Total_Cost) '</strong>']);
        return
    end
    
    %Find the cost of each transition from observed states as the initial state (the states who are stable)
    [Gate_Cost_Table,Pr_Input_State_Mtx,Next_Initial_State,Next_Dest_State,Cost] = Find_Minimum_Cost_Transition(Pr_Input_State_Mtx,W_D,W_O,SPSD(step).SPSD_Observed_States_Numbers,SPSD_Stable_State_Numbers,Num_var,Pr_Gate_Mtx);
    Total_Cost = Total_Cost+Cost;
    
    display('-------------------------------')
    fprintf(['<strong> ' num2str(Next_Initial_State) '   ->   ' num2str(Next_Dest_State) '</strong>\n'])
    display('-------------------------------')
    display(Gate_Cost_Table)
    
    
end
