clc,clear,close all

%-----load parameters-----
load('human_parameter.mat')

vaccine_varient=[1,1,2,2];
vaccine_type=[2,2,2,2];
vaccine_amount=[30,30,30,30];
vaccine_time_pre=[10,40];

distance=(0:0.1:8);
self_neu=ones(2,1);

varient_color=[204,51,17;51,186,237]/255;

vaccine_correspond=[1;2];
gap_sample=14;

interval=(1:36)*30;

Ab_1_sample=zeros(length(interval),length(distance));
Ab_2_sample=zeros(length(interval),length(distance));
M_1_sample=zeros(length(interval),length(distance));
M_2_sample=zeros(length(interval),length(distance));
titer_1_sample=zeros(length(interval),length(distance));
titer_2_sample=zeros(length(interval),length(distance));

Ag_interval_dis=cell(length(interval),length(distance));
Ab_interval_dis=cell(length(interval),length(distance));
F_interval_dis=cell(length(interval),length(distance));
M_off_interval_dis=cell(length(interval),length(distance));
M_on_interval_dis=cell(length(interval),length(distance));
t_interval_dis=cell(length(interval),length(distance));
L=length(distance);


p = gcp('nocreate');
if isempty(p)
    disp('启动并行运算，核心数：12');
    parpool('local', 12);
else
    disp(['并行运算已启动，核心数：' num2str(p.NumWorkers)]);
end

parfor i=(1:length(interval))
     vaccine_time=[vaccine_time_pre,vaccine_time_pre+30+interval(i)]
     i
     for j=(1:L)
        X=eye(2,2);
        cross_neu=2^(-distance(j));
        X(1,2)=cross_neu;
        X(2,1)=cross_neu;      
        [t,Rt,Agt,Abt,Ft,M_offt,M_ont]=simIPD_last(vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,X,self_neu,parameter_set,gap_sample);
        Mt=M_offt+M_ont;
        titert_2=X(2,:)'.*Abt.*self_neu;
        titert_1=X(1,:)'.*Abt.*self_neu;
        t_sample=vaccine_time(end)+gap_sample;
        [Min,ind]=min(abs(t-t_sample'),[],2);
        titer_2_sample(i,j)=sum(titert_2(:,ind));
        titer_1_sample(i,j)=sum(titert_1(:,ind));
        Ab_1_sample(i,j)=Abt(1,ind);
        Ab_2_sample(i,j)=Abt(2,ind);
        M_1_sample(i,j)=Mt(1,ind);
        M_2_sample(i,j)=Mt(2,ind);     
        Ag_interval_dis{i,j}=Agt;
        Ab_interval_dis{i,j}=Abt;
        F_interval_dis{i,j}=Ft;
        M_off_interval_dis{i,j}=M_offt;
        M_on_interval_dis{i,j}=M_ont;
        t_interval_dis{i,j}=t;
     end
end
%%
save('fig6b_data.mat','Ab_1_sample','Ab_2_sample','M_1_sample','M_2_sample','titer_1_sample','titer_2_sample','Ag_interval_dis','Ab_interval_dis','F_interval_dis','M_off_interval_dis','M_on_interval_dis','t_interval_dis')
