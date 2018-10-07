function [Boolean_Function_Common_1st,SPSD_Observed_States_Numbers] = Generate_Boolean_Function_Common_1st_Step(SPSD_Observed_States_Numbers,Output_States,Num_State,Removed_Function_Index)
%This function finds the Boolean functions having a common set of transitions

%Find the state number of the observed and unobserved states
SPSD_UnObserved_States_Numbers = setdiff(0:Num_State-1,SPSD_Observed_States_Numbers);

SPSD_Observed_States_Numbers = sort(SPSD_Observed_States_Numbers);

%generate the truth tables of the Boolean functions
Boolean_Function_Common_1st= zeros(Num_State,2^length(SPSD_UnObserved_States_Numbers));
Boolean_Function_Common_1st(SPSD_UnObserved_States_Numbers+1,:)= de2bi(0:2^(Num_State-length(SPSD_Observed_States_Numbers))-1,'left-msb' )';
Boolean_Function_Common_1st(SPSD_Observed_States_Numbers+1,:)=repmat(Output_States,1,size(Boolean_Function_Common_1st,2));

%Remove Boolean function which should not be included in the search space
Boolean_Func_Index = bi2de(Boolean_Function_Common_1st','left-msb');

Removed_indx = ismember(Boolean_Func_Index,Removed_Function_Index);
Boolean_Function_Common_1st = Boolean_Function_Common_1st(:,~Removed_indx);

end