function order_index=order_vaccination(batch_vaccine_varient,batch_vaccine_time,batch_vaccine_type,batch_vaccine_amount,titer_group,titer_sample)

%根据最后接种的varient排序
n_titer=length(titer_group);
for i=(1:n_titer)
    group=titer_group(i);
    sample=titer_sample(i);
    
    vaccine_varient_i= batch_vaccine_varient(group,1:sample+1);
    vaccine_time_i= batch_vaccine_time(group,1:sample+1);
    vaccine_type_i= batch_vaccine_type(group,1:sample+1);
    vaccine_amount_i= batch_vaccine_amount(group,1:sample+1);
    
    first_type(i)=vaccine_type_i(1);
    end_varient(i)=vaccine_varient_i(end);
    end_type(i)=vaccine_type_i(end);
    injection_times(i)=length(vaccine_varient_i);
    vaccine_time_sum(i)=sum(vaccine_time_i);
    vaccine_amount_sum(i)=sum(vaccine_amount_i);
end
weight=[100000000000,100000000,1000000,10000,100,0.1];
order_inform=weight(1)*first_type+weight(2)*end_varient+weight(3)*end_type+weight(4)*injection_times+weight(5)*vaccine_amount_sum+weight(6)*vaccine_time_sum;

[order,order_index]=sort(order_inform);




end
