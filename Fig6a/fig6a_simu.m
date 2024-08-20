clc,clear,close all

%-----load parameters-----
load('human_parameter.mat')

vaccine_varient=[1,1,2,2];
vaccine_type=[2,2,2,2];
vaccine_time=[10,40,400,430];
booster_dosage=30;

X=eye(2,2);
self_neu=ones(2,1);

pre_varient=1;
gap_sample=14;

pre_amount=(1:50);
distance=(0:0.1:10); %antigenic distance
% pre_amount=(1:2);
% distance=(0:0.1:0.2); %antigenic distance

%-----define cell for saving---------
Ag_amount_dis=cell(length(pre_amount),length(distance));
Ab_amount_dis=cell(length(pre_amount),length(distance));
F_amount_dis=cell(length(pre_amount),length(distance));
M_off_amount_dis=cell(length(pre_amount),length(distance));
M_on_amount_dis=cell(length(pre_amount),length(distance));
t_amount_dis=cell(length(pre_amount),length(distance));


p = gcp('nocreate');
if isempty(p)
    disp('启动并行运算，核心数：12');
    parpool('local', 12);
else
    disp(['并行运算已启动，核心数：' num2str(p.NumWorkers)]);
end
pre_vac=2;

pre_vaccine_varient=vaccine_varient(1:pre_vac);
pre_vaccine_type=vaccine_type(1:pre_vac);
pre_vaccine_time=vaccine_time(1:pre_vac);
gap=vaccine_time(pre_vac+1)-vaccine_time(pre_vac);
gap
pre_X=X(1,1);
pre_self_neu=self_neu(1);


bo_vaccine_varient=vaccine_varient(pre_vac+1:end);
bo_vaccine_type=vaccine_type(pre_vac+1:end);
bo_vaccine_time=vaccine_time(pre_vac+1:end);

L=length(distance);

%%
tic
parfor i=(1:length(pre_amount))
    i
    vaccine_amount=[pre_amount(i),pre_amount(i),booster_dosage,booster_dosage];
    %-------background_simulation-------

    pre_vaccine_amount=vaccine_amount(1:pre_vac);

    [pre_t,pre_Rt,pre_Agt,pre_Abt,pre_Ft,pre_M_offt,pre_M_ont]=simIPD_last(pre_vaccine_varient,pre_vaccine_time,pre_vaccine_type,pre_vaccine_amount,pre_X,pre_self_neu,parameter_set,gap);

    bo_vaccine_amount=vaccine_amount(pre_vac+1:end);

    for j=(1:L)
        para_X=eye(2,2);
        cross_neu=2^(-distance(j));
        para_X(1,2)=cross_neu;
        para_X(2,1)=cross_neu;   
        [t,Rt,Agt,Abt,Ft,M_offt,M_ont]=simIPD_continue(bo_vaccine_varient,bo_vaccine_time,bo_vaccine_type,bo_vaccine_amount,para_X,self_neu,parameter_set,gap_sample,pre_t,pre_Rt,pre_Agt,pre_Abt,pre_Ft,pre_M_offt,pre_M_ont,pre_varient);
        Ag_amount_dis{i,j}=Agt;
        Ab_amount_dis{i,j}=Abt;
        F_amount_dis{i,j}=Ft;
        M_off_amount_dis{i,j}=M_offt;
        M_on_amount_dis{i,j}=M_ont;
        t_amount_dis{i,j}=t;
        
    end
end
toc
save('fig6a_data.mat','pre_amount','distance','t_amount_dis','Ag_amount_dis','Ab_amount_dis','F_amount_dis','M_off_amount_dis','M_on_amount_dis')