
clc
clear SPSD;
load('Removed_Function_Index.mat');
%% Find the bit combinations (states) of input bits
%input simulation

%% ABC
SPSD(1).Initial_State_str = '011';
SPSD(1).Target_State_str = '0';
SPSD(1).Perturbed_Input_index_str = '1 2 3';
SPSD(1).Output_States_After_Pert_str = '100';

SPSD(2).Initial_State_str = '111';
SPSD(2).Perturbed_Input_index_str = '2 3';
SPSD(2).Output_States_After_Pert_str = '00';

SPSD(3).Initial_State_str = '101';
SPSD(3).Perturbed_Input_index_str = '2 3';
SPSD(3).Output_States_After_Pert_str = '00';

SPSD(4).Initial_State_str = '100';
SPSD(4).Perturbed_Input_index_str = '1';
SPSD(4).Output_States_After_Pert_str = '0';

%This variable shows which input state combinations are stable (could be used as the initial state for transitions)
SPSD_Stable_State_Numbers = [1 1 1 1 1 1 1 1];


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
% SPSD(1).Target_State_str = '1';
% SPSD(1).Perturbed_Input_index_str = '1 2 3 4';
% SPSD(1).Output_States_After_Pert_str = '0000';
% 
% 
% SPSD(2).Initial_State_str = '1110';
% SPSD(2).Perturbed_Input_index_str = '1 2 3';
% SPSD(2).Output_States_After_Pert_str = '000';
% 
% SPSD(3).Initial_State_str = '1100';
% SPSD(3).Perturbed_Input_index_str = '1 2 4';
% SPSD(3).Output_States_After_Pert_str = '000';
% 
% SPSD(4).Initial_State_str = '0100';
% SPSD(4).Perturbed_Input_index_str = '1 2 3 4';
% SPSD(4).Output_States_After_Pert_str = '0000';
% 
% SPSD(5).Initial_State_str = '0000';
% SPSD(5).Perturbed_Input_index_str = '1 2 3 4';
% SPSD(5).Output_States_After_Pert_str = '0000';
% 
% SPSD(6).Initial_State_str = '0010';
% SPSD(6).Perturbed_Input_index_str = '4';
% SPSD(6).Output_States_After_Pert_str = '0';
% 
% SPSD(7).Initial_State_str = '1011';
% SPSD(7).Perturbed_Input_index_str = '3';
% SPSD(7).Output_States_After_Pert_str = '0';

%% Initialization
Num_var = length(SPSD(1).Initial_State_str) ;
Num_State = 2^Num_var;

Gate(2).gate = {'~A~B'; '~AB';'A~B';'AB'};
Gate(3).gate = {'~A~B~C';'~A~BC';'~AB~C';'~ABC';'A~B~C';'A~BC';'AB~C';'ABC'};
Gate(4).gate = {'~A~B~C~D';'~A~B~CD';'~A~BC~D';'~A~BCD';'~AB~C~D';'~AB~CD';'~ABC~D';'~ABCD';'A~B~C~D';'A~B~CD';'A~BC~D';'A~BCD';'AB~C~D';'AB~CD';'ABC~D';'ABCD'};

%% Get Data from user

% % %Get the steady state of the TFs
% % Initial_State_str = input('Enter the steady state of TFs in binary:','s');
% % if(length(Initial_State_str) ~= Num_var)
% %     disp('The length of the input for steady state is WRONG!');
% %     return
% % end
SPSD(1).Initial_State_Binary = SPSD(1).Initial_State_str=='1';

% % %Get the steady state of the target gene
% % Target_State_str = input('Enter the steady state of the target gene in binary:','s');
% % if(length(Target_State_str) ~= 1)
% %     disp('The length of the input should be 1!');
% %     return
% % end
Target_State_Binary_1 = SPSD(1).Target_State_str=='1';

% % %Get the index of the perturbed inputs in the first step
% % SPSD(1).Perturbed_Input_index_str = input('Enter the index of perturbed inpute in order with space(e.g., 1 2 3 or 1 3):','s');
SPSD(1).Perturbed_Input_index = str2num(SPSD(1).Perturbed_Input_index_str);

% % if(length(SPSD(1).Perturbed_Input_index) > Num_var || isempty(SPSD(1).Perturbed_Input_index) )
% %     disp(['The length of  input should not be between 1 and ' num2str(Num_var)]);
% %     return
% % end
% % 
% % if(max(SPSD(1).Perturbed_Input_index)>Num_var)
% %     disp(['values should be less than ' num2str(Num_var) '.'] );
% %     return
% % end



%Get the transition in the target output by doing knockdown and overexperission
% % SPSD(1).Output_States_After_Pert_str = [];
% % while (length(SPSD(1).Output_States_After_Pert_str) ~= length(SPSD(1).Perturbed_Input_index))
% %     SPSD(1).Output_States_After_Pert_str = '000';%input('Enter the perturbed outputs:','s');
% %     if(length(SPSD(1).Output_States_After_Pert_str) ~= length(SPSD(1).Perturbed_Input_index))
% %         disp(['The length of the input for output states should be ' num2str(length(SPSD(1).Perturbed_Input_index))]);
% %     end
% % end
SPSD(1).Output_States_After_Pert = SPSD(1).Output_States_After_Pert_str=='1';
SPSD(1).Output_States = [Target_State_Binary_1;SPSD(1).Output_States_After_Pert'];


%Find the new input states after perturbation (the first row of SPSD_States is
%the initial steady state of the inputs)
SPSD(1).SPSD_Perturbed_States = repmat(SPSD(1).Initial_State_Binary,length(SPSD(1).Perturbed_Input_index)+1,1);
for i = 1:length(SPSD(1).Perturbed_Input_index)
    indx = SPSD(1).Perturbed_Input_index(i);
    SPSD(1).SPSD_Perturbed_States(i+1,indx) = ~SPSD(1).SPSD_Perturbed_States(i+1,indx);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                First Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate the truth table of Boolean function having common transitions %(1st step)
[SPSD(1).Boolean_Function_Common , SPSD_Observed_States_Numbers] = Generate_Boolean_Function_Common_1(SPSD(1).SPSD_Perturbed_States,SPSD(1).Output_States,Num_State,SPSD(1).Output_States_After_Pert,Target_State_Binary_1,Removed_Function_Index(Num_var).index);

%Find the probability of each gate by counting the number of the ones in truth table for each input state
[Pr_Input_State_Mtx , Pr_Gate_Mtx] = Find_Gate_Probability(SPSD(1).Boolean_Function_Common,Gate,Num_var,1);

%Set the cost weight of knock-down (D) and over-expression (O)
W_D = 1;
W_O = 2;

Statble_Input_State = SPSD(1).Initial_State_Binary;
[T,Pr_Input_State_Mtx] = Find_Minimum_Cost_Transition(Pr_Input_State_Mtx,W_D,W_O,SPSD_Observed_States_Numbers,SPSD_Stable_State_Numbers,Num_var,Pr_Gate_Mtx);

display(T)


return
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                2nd Step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Get the 2nd initial steady state 
% % SPSD(2).Initial_State_str = input('Enter the 2nd initial steady state of TFs in binary:','s');
% % if(length(SPSD(2).Initial_State_str) ~= Num_var)
% %     disp('The length of the input for steady state is WRONG!');
% %     return
% % end
SPSD(2).Initial_State_Binary = SPSD(2).Initial_State_str=='1';


% % %Get the index of the perturbed inputs in the 2nd step
% % SPSD(2).Perturbed_Input_index_str = input('Enter the index of perturbed inpute in order with space(e.g., 1 2 3 or 1 3):','s');
SPSD(2).Perturbed_Input_index = str2num(SPSD(2).Perturbed_Input_index_str);

% % if(length(SPSD(2).Perturbed_Input_index) > Num_var || isempty(SPSD(2).Perturbed_Input_index) )
% %     disp(['The length of  input should not be between 1 and ' num2str(Num_var)]);
% %     return
% % end
% % 
% % if(max(SPSD(2).Perturbed_Input_index)>Num_var)
% %     disp(['values should be less than ' num2str(Num_var) '.'] );
% %     return
% % end

% % Get the transition in the target output by doing knockdown and overexperission
% % SPSD(2).Output_States_After_Pert_str = [];
% % while (length(SPSD(2).Output_States_After_Pert_str) ~= length(SPSD(2).Perturbed_Input_index))
% %     SPSD(2).Output_States_After_Pert_str = '000';%input('Enter the perturbed outputs:','s');
% %     if(length(SPSD(2).Output_States_After_Pert_str) ~= length(SPSD(2).Perturbed_Input_index))
% %         disp(['The length of the input for output states should be ' num2str(length(SPSD(2).Perturbed_Input_index))]);
% %     end
% % % % end

for step = 2:length(SPSD)

    SPSD(step).Initial_State_Binary = SPSD(step).Initial_State_str=='1';
    SPSD(step).Perturbed_Input_index = str2num(SPSD(step).Perturbed_Input_index_str);
    SPSD(step).Output_States_After_Pert = SPSD(step).Output_States_After_Pert_str=='1';
    
    
    %Generate the truth table of Boolean function having common transitions %(2nd step)
    SPSD = Generate_Boolean_Function_Common_2(SPSD,step);
    
    %Find the probability of each gate by counting the number of the ones in truth table for each input state
    Pr_Input_State_Mtx = Find_Gate_Probability(SPSD(step).Boolean_Function_Common,Gate,Num_var,step);

end
