function Sub_State = Creat_Sub_State(States,Visible_TF_Effect_Index,Data_Type,index,state_num)

    switch Data_Type
        case 'D'
            Sub_State = States(index,Visible_TF_Effect_Index.D(state_num).Output_index);
        case 'DO'
            Sub_State = States(index,Visible_TF_Effect_Index.DO(state_num).Output_index);
        case 'D_PPI'
            Sub_State = States(index,Visible_TF_Effect_Index.D_PPI(state_num).Output_index);
        case 'DO_PPI'
            Sub_State = States(index,Visible_TF_Effect_Index.DO_PPI(state_num).Output_index);
        case {'D_cis','DO_cis'}
            Sub_State = States(index,:);
    end
end