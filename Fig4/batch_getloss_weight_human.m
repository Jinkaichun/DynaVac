function [loss,n_inverse]=batch_getloss_weight_human(batch_titer_lgFC,batch_titer_est_lgFC,weight)

%term_weight=ones(size(batch_titer_lgFC)); %每一项的权重

%loss1=sum(sum(term_weight.*(batch_titer_est_lgFC-batch_titer_lgFC).^2));

mean_lgFC=mean(batch_titer_lgFC(:),'omitnan');
TSS=sum((batch_titer_lgFC(:)-mean_lgFC).^2,'omitnan');
ESS=sum((batch_titer_lgFC(:)-batch_titer_est_lgFC(:)).^2,'omitnan');

loss1=ESS/TSS;



n_varient=size(batch_titer_lgFC,1);

for i=(1:n_varient)
    batch_titer_est_vect=batch_titer_est_lgFC(i,:);
    [batch_titer_vect_sort,sort_ind]=sort(batch_titer_lgFC(i,:));
    batch_titer_est_sort=batch_titer_est_vect(sort_ind);
    n(i)=inverse_num(batch_titer_est_sort);
    
end

%loss2=sum(n); %将针对每一种varient的逆序数相加，获得loss


n_titer=size(batch_titer_lgFC,2);
loss2=mean(n/(n_titer*(n_titer-1)/2));

batch_titer_est_vect=batch_titer_est_lgFC(:);
[batch_titer_vect_sort,sort_ind]=sort(batch_titer_lgFC(:));
batch_titer_est_sort=batch_titer_est_vect(sort_ind);

n_all=inverse_num(batch_titer_est_sort);
n_titer_all=length(batch_titer_est_sort);

loss3=n_all/(n_titer_all*(n_titer_all-1)/2);



loss=weight(1)*loss1+weight(2)*loss2+weight(3)*loss3;
n_inverse=n;
R=1-loss1
rank=(1-n/(n_titer*(n_titer-1)/2))*100
rank_all=(1-n_all/(n_titer_all*(n_titer_all-1)/2))*100
%[loss1,loss2,loss3,loss];%display


end
