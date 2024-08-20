function [batch_t,batch_Rt,batch_Agt,batch_Abt,batch_Ft,batch_Mofft,batch_Mont]=batch_simIPD(batch_vaccine_varient,batch_vaccine_time,batch_vaccine_type,batch_vaccine_amount,X,self_neu,parameter_set,vaccine_correspond)

batch_num=size(batch_vaccine_varient,1);
batch_t=cell(batch_num,1);
batch_Rt=cell(batch_num,1);
batch_Agt=cell(batch_num,1);
batch_Abt=cell(batch_num,1);
batch_Ft=cell(batch_num,1);
batch_Mofft=cell(batch_num,1);
batch_Mont=cell(batch_num,1);


for batch=(1:batch_num)
    
    vaccine_varient=batch_vaccine_varient(batch,:);
    vaccine_time=batch_vaccine_time(batch,:);
    vaccine_type=batch_vaccine_type(batch,:);
    vaccine_amount=batch_vaccine_amount(batch,:);

    eff_ind=vaccine_varient~=0;
    vaccine_varient= vaccine_varient(eff_ind);
    vaccine_time= vaccine_time(eff_ind);
    vaccine_type= vaccine_type(eff_ind);
    vaccine_amount= vaccine_amount(eff_ind);
    
    vaccine_varient=vaccine_correspond(vaccine_varient,:)';
    
    
    sub_varient=unique(vaccine_varient);
    sub_varient=sub_varient(sub_varient~=0);
    sub_X=X(sub_varient,sub_varient);
    sub_self_neu=self_neu(sub_varient);
    
    [t,Rt,Agt,Abt,Ft,Mofft,Mont]=simIPD(vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,sub_X,sub_self_neu,parameter_set);
    
    batch_t{batch}=t;
    batch_Rt{batch}=Rt;
    batch_Agt{batch}=Agt;
    batch_Abt{batch}=Abt;
    batch_Ft{batch}=Ft;
    batch_Mofft{batch}=Mofft;
    batch_Mont{batch}=Mont;

end