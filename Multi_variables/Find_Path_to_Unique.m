clc
load('Combination_Family.mat');

Num_var = 3 ;
Num_State = 2^Num_var;

%% Find the bit combinations (states) of input bits
Input_States = de2bi(Num_State-1:-1:0,'left-msb');

%% Find the hamming distances between different states

%The (i,j)th element shows the hamming distance between states number i & j
Hamming_Dist = zeros(Num_State);
for i=1:Num_State-1
    for j=i+1:Num_State
        Hamming_Dist (i,j) = pdist2(Input_States(i,:),Input_States(j,:),'Hamming');
    end
end
Hamming_Dist = (Hamming_Dist+Hamming_Dist')*Num_var;

%% Get Data from user

%Get the steady state of the TFs
Steady_State_str = input('Enter the steady state of TFs in binary:','s');
if(length(Steady_State_str) ~= Num_var)
    disp('The length of the input for steady state is WRONG!');
    return
end
Steady_State_Binary = Steady_State_str=='1';
Steady_State_Number = Num_State-bi2de(Steady_State_Binary,'left-msb');

%Get the steady state of the target gene
Target_State_str = input('Enter the steady state of the target gene in binary:','s');
if(length(Target_State_str) ~= 1)
    disp('The length of the input should be 1!');
    return
end
Target_State_Binary = Target_State_str=='1';


%Get the transition in the target output by doing knockdown and overexperission
State_Transition_str = input('Enter the target transitions (D+O) in binary:','s');
if(length(State_Transition_str) ~= Num_var)
    disp('The length of the input for transitions is WRONG!');
    return
end
State_Transition_Binary = State_Transition_str=='1';

Transition = [Target_State_Binary State_Transition_Binary]';


Neighbor_State_Index =  [Steady_State_Number find(CF(Num_var).Interact(Steady_State_Number,:))];
Sub_State = CF(Num_var).States(Neighbor_State_Index,:);

Sate_Index = find(sum(abs(Sub_State-repmat(Transition,1,size(Sub_State,2))))==0);
if length(Sate_Index)==1
    %We are in the unique state
    disp('The gate combination family is as follows:')
    disp(CF(Num_var).States(Num_State+(1:((2^(Num_var-2)))),Sate_Index));
    return
end

%% Find the most probable states to be unique
All_Possible_Combination_Index = CF(Num_var).Gate_Combination_Unique_Transition_indx(:,Sate_Index);
All_Possible_Combination_Index = sort(All_Possible_Combination_Index(:))';

[freq,state_num]=hist(All_Possible_Combination_Index,1:Num_State);

[Sorted_State_Freq, All_Possible_Combination_Index]=sort(freq,'descend');
Sorted_State_Freq = Sorted_State_Freq(Sorted_State_Freq~=0);
All_Possible_Combination_Index = All_Possible_Combination_Index(Sorted_State_Freq~=0);

for i=max(Sorted_State_Freq):-1:1
   freq_ind = find( Sorted_State_Freq == i);
   if isempty(freq_ind)
       continue;
   end
   State_ind = All_Possible_Combination_Index(freq_ind);
   Dist = Hamming_Dist(Steady_State_Number,State_ind);
   [~ , Dist_Sort_Indx] = sort(Dist);
   
   All_Possible_Combination_Index(freq_ind) = All_Possible_Combination_Index(freq_ind(Dist_Sort_Indx));
end
All_Possible_Combination_Index = [All_Possible_Combination_Index ;Hamming_Dist(Steady_State_Number,All_Possible_Combination_Index);round(Sorted_State_Freq/sum(Sorted_State_Freq)*100)]

 Row_Name = {'Unique state';'Hamming Dist.';'Probability'};
