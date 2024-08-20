function n=inverse_num(A)
[inda,indb]=find(triu(ones(length(A)),1)); %triu-取上三角矩阵
n=sum(((A(inda)-A(indb))>0));
end