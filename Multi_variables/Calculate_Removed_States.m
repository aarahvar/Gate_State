function Removed_States = Calculate_Removed_States(Num_var,Target_Output,Removed_States,Indexes)
for num_var = 2:Num_var-1
    for i=1:size(Target_Output(num_var).States,2)
        g= Target_Output(num_var).States(:,i)';
        g=repmat(g,2^(Num_var-num_var),1);
        g=g(:);
        Removed_States = [Removed_States g(Indexes)];
    end
end

end

