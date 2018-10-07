function g_index = Get_Gate_Index(G_str,Num_Gate)
k=0;
for i=1:length(G_str)
    if G_str(i)~='.'
        k=k+1;
        G_Num = str2num(G_str(i));
        if(G_Num)
            g_index(k) = G_Num;
        else
            g_index(k) = Num_Gate;
        end
        
    end
end
end