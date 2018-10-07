function Num_State_Members = Get_Num_State_Members(Visible_TF_Effect_Index,Data_Type,Num_State)

Num_State_Members = zeros(1,Num_State);
for k=1:Num_State
    switch Data_Type
        case 'D'
            Num_State_Members(k) = length(Visible_TF_Effect_Index.D(k).Output_index);
        case 'DO'
            Num_State_Members(k) = length(Visible_TF_Effect_Index.DO(k).Output_index);
        case 'D_PPI'
            Num_State_Members(k) = length(Visible_TF_Effect_Index.D_PPI(k).Output_index);
        case 'DO_PPI'
            Num_State_Members(k) = length(Visible_TF_Effect_Index.DO_PPI(k).Output_index);
        case {'D_cis','DO_cis'}
            Num_State_Members(k) = length(Visible_TF_Effect_Index.cis(k).Output_index);
    end
end

end