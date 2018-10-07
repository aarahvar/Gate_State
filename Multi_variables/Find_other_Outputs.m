Input_States = de2bi(0:2^(2^Num_var)-1,'left-msb')';

ABC = States(1:8,:);
BCA = States([1,3,5,7,2,4,6,8],:);
ACB = States([1,3,2,4,5,7,6,8],:);

Combinations = [ABC BCA ACB];

Others = [];
for i=1:size(Input_States,2)
   indx = find(sum(abs(Combinations-repmat(Input_States(:,i),1,size(Combinations,2))))==0); 
   
   if isempty(indx)
       Others = [Others Input_States(:,i)];
   end
    
end