function [Pr_Input_State_Mtx , Pr_Gate_Mtx] = Find_Gate_Probability(Boolean_Function_Common,Gate,Num_var,Step)
%This functions calculates the probability of occurance of each gate combintion by counting the probability of ocurance of 1s in the truth tables
%The output is a matrix that the first column is the list of sorted probabilities and other columns in each row is the binary state of the corresponding gate inpute

%Calculate the probability of occurance of 1s for each gate combination
Pr_Gate = sum(Boolean_Function_Common,2)/size(Boolean_Function_Common,2);

%sort the probabilities
[Pr_Gate_Sort,Pr_Gate_Sort_Index] = sort(Pr_Gate,'descend');


%generate the probability-gate string and matrix
Pr_Input_State_Mtx = [];
Pr_Gate_Mtx = [];
Pr_Gate_str = 'Pr(A,B,C)=';
for pr=find(Pr_Gate_Sort)'
    g=Pr_Gate_Sort_Index(pr);
    if Pr_Gate(g)==1
        if Pr_Gate_str(end)=='='
            Pr_Gate_str = [Pr_Gate_str ' ' cell2mat(Gate(Num_var).gate(g))];
        else
            Pr_Gate_str = [Pr_Gate_str ' + ' cell2mat(Gate(Num_var).gate(g))];
        end
    else
        if Pr_Gate_str(end)=='='
            Pr_Gate_str = [Pr_Gate_str ' (' num2str(Pr_Gate(g),2) ')' cell2mat(Gate(Num_var).gate(g))];
        else
            Pr_Gate_str = [Pr_Gate_str ' + (' num2str(Pr_Gate(g),2) ')' cell2mat(Gate(Num_var).gate(g))];
        end
        
    end
    
    Pr_Input_State_Mtx = [Pr_Input_State_Mtx;Pr_Gate_Sort(pr) dec2bin(g-1,Num_var)=='1'];
    Pr_Gate_Mtx{end+1} = (Gate(Num_var).gate(g));
    
end


%Print the probability-gate output
% disp('**************************************************************************************');

Pr_str = ['<strong>Step '  num2str(Step-1) '</strong>:    ' Pr_Gate_str];
star = [];
for i=1:length(Pr_str)-15
    star = [star '*'];
end

disp(star);
fprintf(Pr_str);
fprintf('\n');
disp(star);



end