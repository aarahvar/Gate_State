function SPSD = Generate_Boolean_Function_Common_2(SPSD,step,Output_State)
%This function finds the Boolean functions having a common set of transitions for after the first step

Initial_State_Binary_2       = SPSD(step).Initial_State_Binary;
Boolean_Function_Common_1   = SPSD(step-1).Boolean_Function_Common;

SPSD_Observed_States_1     = SPSD(step-1).SPSD_Observed_States;
Output_States_1             = SPSD(step-1).Output_States;

SPSD_Observed_States_2 = [SPSD_Observed_States_1;Initial_State_Binary_2];
Output_States_2 = [Output_States_1;Output_State(bi2de(Initial_State_Binary_2,'left-msb')+1)];


%Find the state number of the observed and unobserved states
SPSD_Observed_States_Numbers_2 = bi2de(SPSD_Observed_States_2,'left-msb');
[SPSD_Observed_States_Numbers_2, unique_indx] = unique(SPSD_Observed_States_Numbers_2,'first');
Output_States_2 = Output_States_2(unique_indx);
SPSD_Observed_States_2 = SPSD_Observed_States_2(unique_indx,:);
% [SPSD_Observed_States_Numbers_2,sort_indx] = sort(SPSD_Observed_States_Numbers_2);
% Output_States_2 = Output_States_2(sort_indx);


Observed_States = Boolean_Function_Common_1(SPSD_Observed_States_Numbers_2+1,:);
Common_index = sum(abs(Observed_States-repmat(Output_States_2,1,size(Boolean_Function_Common_1,2))))==0;
Boolean_Function_Common_2 = Boolean_Function_Common_1(:,Common_index);

SPSD(step).Boolean_Function_Common = Boolean_Function_Common_2;
SPSD(step).SPSD_Observed_States = SPSD_Observed_States_2;
SPSD(step).Output_States = Output_States_2;
SPSD(step).SPSD_Observed_States_Numbers = SPSD_Observed_States_Numbers_2;
end