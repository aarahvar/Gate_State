clc
Struct.Mean_Common_State = [];
Struct.Num_State_Members = [];
Struct.index = [];
Struct.Visibility_Unambiguity_Score = [];
Struct.Mean_Common_in_Visible = [];
load('Target_Output_Full.mat');

SNames = fieldnames(Common_State_Visible);
for Num_var=2:4
    a=[];

    States = Target_Output(Num_var).States;
    for loopIndex = 1:numel(SNames)
        if isempty(Common_State_Visible(Num_var).(SNames{loopIndex}))
            continue;
        end
        Struct = Common_State_Visible(Num_var).(SNames{loopIndex});
        P_Visibility_in_State = Struct.Num_State_Members/size(States,2);
        ind = (Struct.Mean_Common_State~=0);
        
        Struct= rmfield(Struct,'Visibility_Unambiguity_Score');
        Struct.Unambiguity_Score = sum(Struct.Mean_Common_State(ind)./P_Visibility_in_State(ind))
        a(end+1)=Struct.Unambiguity_Score;
    end
end
plot(a)
% for Num_var=4:-1:2
% Common_State_Visible(Num_var).Single_D.Mean_Common_State = [];
% Common_State_Visible(Num_var).Single_D.Num_State_Members = [];
% Common_State_Visible(Num_var).Single_D.index = [];
% Common_State_Visible(Num_var).Single_D.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Single_D.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Single_D_PPI.Mean_Common_State = [];
% Common_State_Visible(Num_var).Single_D_PPI.Num_State_Members = [];
% Common_State_Visible(Num_var).Single_D_PPI.index = [];
% Common_State_Visible(Num_var).Single_D_PPI.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Single_D_PPI.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Single_D_cis.Mean_Common_State = [];
% Common_State_Visible(Num_var).Single_D_cis.Num_State_Members = [];
% Common_State_Visible(Num_var).Single_D_cis.index = [];
% Common_State_Visible(Num_var).Single_D_cis.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Single_D_cis.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Double_D.Mean_Common_State = [];
% Common_State_Visible(Num_var).Double_D.Num_State_Members = [];
% Common_State_Visible(Num_var).Double_D.index = [];
% Common_State_Visible(Num_var).Double_D.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Double_D.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Double_D_PPI.Mean_Common_State = [];
% Common_State_Visible(Num_var).Double_D_PPI.Num_State_Members = [];
% Common_State_Visible(Num_var).Double_D_PPI.index = [];
% Common_State_Visible(Num_var).Double_D_PPI.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Double_D_PPI.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Double_D_cis.Mean_Common_State = [];
% Common_State_Visible(Num_var).Double_D_cis.Num_State_Members = [];
% Common_State_Visible(Num_var).Double_D_cis.index = [];
% Common_State_Visible(Num_var).Double_D_cis.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Double_D_cis.Mean_Common_in_Visible = [];
% 
% 
% if Num_var~=2
% Common_State_Visible(Num_var).Triple_D.Mean_Common_State = [];
% Common_State_Visible(Num_var).Triple_D.Num_State_Members = [];
% Common_State_Visible(Num_var).Triple_D.index = [];
% Common_State_Visible(Num_var).Triple_D.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Triple_D.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Triple_D_PPI.Mean_Common_State = [];
% Common_State_Visible(Num_var).Triple_D_PPI.Num_State_Members = [];
% Common_State_Visible(Num_var).Triple_D_PPI.index = [];
% Common_State_Visible(Num_var).Triple_D_PPI.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Triple_D_PPI.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Triple_D_cis.Mean_Common_State = [];
% Common_State_Visible(Num_var).Triple_D_cis.Num_State_Members = [];
% Common_State_Visible(Num_var).Triple_D_cis.index = [];
% Common_State_Visible(Num_var).Triple_D_cis.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Triple_D_cis.Mean_Common_in_Visible = [];
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Common_State_Visible(Num_var).Single_DO.Mean_Common_State = [];
% Common_State_Visible(Num_var).Single_DO.Num_State_Members = [];
% Common_State_Visible(Num_var).Single_DO.index = [];
% Common_State_Visible(Num_var).Single_DO.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Single_DO.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Single_DO_PPI.Mean_Common_State = [];
% Common_State_Visible(Num_var).Single_DO_PPI.Num_State_Members = [];
% Common_State_Visible(Num_var).Single_DO_PPI.index = [];
% Common_State_Visible(Num_var).Single_DO_PPI.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Single_DO_PPI.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Single_DO_cis.Mean_Common_State = [];
% Common_State_Visible(Num_var).Single_DO_cis.Num_State_Members = [];
% Common_State_Visible(Num_var).Single_DO_cis.index = [];
% Common_State_Visible(Num_var).Single_DO_cis.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Single_DO_cis.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Double_DO.Mean_Common_State = [];
% Common_State_Visible(Num_var).Double_DO.Num_State_Members = [];
% Common_State_Visible(Num_var).Double_DO.index = [];
% Common_State_Visible(Num_var).Double_DO.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Double_DO.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Double_DO_PPI.Mean_Common_State = [];
% Common_State_Visible(Num_var).Double_DO_PPI.Num_State_Members = [];
% Common_State_Visible(Num_var).Double_DO_PPI.index = [];
% Common_State_Visible(Num_var).Double_DO_PPI.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Double_DO_PPI.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Double_DO_cis.Mean_Common_State = [];
% Common_State_Visible(Num_var).Double_DO_cis.Num_State_Members = [];
% Common_State_Visible(Num_var).Double_DO_cis.index = [];
% Common_State_Visible(Num_var).Double_DO_cis.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Double_DO_cis.Mean_Common_in_Visible = [];
% 
% if Num_var~=2
% Common_State_Visible(Num_var).Triple_DO.Mean_Common_State = [];
% Common_State_Visible(Num_var).Triple_DO.Num_State_Members = [];
% Common_State_Visible(Num_var).Triple_DO.index = [];
% Common_State_Visible(Num_var).Triple_DO.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Triple_DO.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Triple_DO_PPI.Mean_Common_State = [];
% Common_State_Visible(Num_var).Triple_DO_PPI.Num_State_Members = [];
% Common_State_Visible(Num_var).Triple_DO_PPI.index = [];
% Common_State_Visible(Num_var).Triple_DO_PPI.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Triple_DO_PPI.Mean_Common_in_Visible = [];
% 
% Common_State_Visible(Num_var).Triple_DO_cis.Mean_Common_State = [];
% Common_State_Visible(Num_var).Triple_DO_cis.Num_State_Members = [];
% Common_State_Visible(Num_var).Triple_DO_cis.index = [];
% Common_State_Visible(Num_var).Triple_DO_cis.Visibility_Unambiguity_Score = [];
% Common_State_Visible(Num_var).Triple_DO_cis.Mean_Common_in_Visible = [];
% end
% end