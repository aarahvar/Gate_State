Num_func = [ 10 218 64594];
load('C:\Users\Amir\MATLAB\Gate_State\Multi_variables\Visible_TF_Effect_Index.mat')
Num_var = 2;
for Num_var=2:4
    for indx = 1:2^Num_var
        Pfvs(Num_var).P_D(indx,:) = [2^Num_var-indx length(Visible_TF_Effect_Index(Num_var).D(indx).Output_index)/Num_func(Num_var-1)];
    end
end

for Num_var=2:4
    for indx = 1:2^Num_var
        Pfvs(Num_var).P_D_PPI(indx,:) = [2^Num_var-indx length(Visible_TF_Effect_Index(Num_var).D_PPI(indx).Output_index)/Num_func(Num_var-1)];
    end
end

for Num_var=2:4
    for indx = 1:2^Num_var
        Pfvs(Num_var).P_DO(indx,:) = [2^Num_var-indx length(Visible_TF_Effect_Index(Num_var).DO(indx).Output_index)/Num_func(Num_var-1)];
    end
end

for Num_var=2:4
    for indx = 1:2^Num_var
        Pfvs(Num_var).P_DO_PPI(indx,:) = [2^Num_var-indx length(Visible_TF_Effect_Index(Num_var).DO_PPI(indx).Output_index)/Num_func(Num_var-1)];
    end
end
