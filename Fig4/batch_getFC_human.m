function [batch_titer_lgFC,batch_titer_est_lgFC,titer_group,titer_time,titer_sample]=batch_getFC_human(batch_vaccine_varient,batch_vaccine_time,batch_t,batch_Abt,batch_sample,titers,X_pse,self_neu,vaccine_correspond)

batch_num=size(batch_vaccine_varient,1);
batch_titer=[];
batch_titer_est=[];
titer_group=[];
titer_time=[];
titer_sample=[];
for batch=(1:batch_num)
    
    t=batch_t{batch};
    t_sample=batch_sample(batch,:);
    t_sample=t_sample(t_sample~=0);
    
    Abt=batch_Abt{batch};

    vaccine_varient=batch_vaccine_varient(batch,:);
    vaccine_time=batch_vaccine_time(batch,:);
    
    eff_ind=vaccine_varient~=0; %提取有效的数据索引
    vaccine_varient= vaccine_varient(eff_ind);

    
    vaccine_time= vaccine_time(eff_ind);
    vaccine_varient=vaccine_correspond(vaccine_varient,:)';
    
    sub_varient=unique(vaccine_varient);
    sub_varient=sub_varient(sub_varient~=0);
    
    sub_X=X_pse(:,sub_varient); %所有行，个别列
    sub_self_neu=self_neu(sub_varient);
    
    titert=sub_X*(Abt.*sub_self_neu);

    [~,ind]=min(abs(t-t_sample'),[],2); %找到时间上最接近的点
    titer_est=titert(:,ind);
    
    titer=titers{batch};

    titer_sample_i=sum(t_sample'>vaccine_time,2)';%每一次采样在第几次接种之后

    
    batch_titer=[batch_titer,titer];
    batch_titer_est=[batch_titer_est,titer_est];
    titer_group=[titer_group,batch*ones(size(t_sample))];%记录每一列数据属于哪个batch
    titer_time=[titer_time,t_sample];
    titer_sample=[titer_sample,titer_sample_i];
    
end

batch_titer_lgFC=log2(batch_titer/batch_titer(1,1));
batch_titer_est_lgFC=log2(batch_titer_est/batch_titer_est(1,1));


end
