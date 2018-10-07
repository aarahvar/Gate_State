clc
clear max_degree;
%Find the states with more changes in the output in state transitions (Change %degree)
for i= 1:size(States,2)
    Change_Degree = sum(Interact.*mod(repmat(States(1:Num_State,i),1,Num_State)+repmat(States(1:Num_State,i),1,Num_State)',2));
    [sorted,index] = sort(Change_Degree,'descend');
    max_degree(i,1) = sorted(1);
    max_degree(i,2) = sorted(2);
end
figure(1);hist(max_degree(:,1))
% 
% States = CF(2).States;
% index = [];
% for i= 1:size(States,2)
%    a=dec2bin((2^Num_var)-Gate_Combination_Unique_Transition_indx(:,i)) ;
%    for m=1:2
%        for n=m+1:3
%            index = [index find(a(m,:)~=a(n,:))];
%            
%        end
%    end
%     
%     
% end
%  figure(1);plot(sort(index))
%  sum(find(index==1)~=0)
%   sum(find(index==2)~=0)