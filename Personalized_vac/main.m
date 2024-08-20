clc,clear,close all
load('X_yunda.mat')
load('self_neu_yunda.mat');
self_neu=self_neu';
load('mix_parameters.mat')
X(5,1)=X(5,1)*sqrt(parameter_set(25));
X(6,1)=X(6,1)*parameter_set(28);
X(6,2)=X(6,2)*parameter_set(30);
X(6,3)=X(6,3)*parameter_set(31);
X(6,4)=X(6,4)*parameter_set(32);

%%
%-----------personalized vaccination strategy---------
vaccine_varient=[1,4,3,2,5;0,0,5,0,6];%V: vaccination variant
%(1: WT, 2: Alpha/Beta, 3: Delta, 4:BA.1, 5: BA.2/BA.5, 6: XBB.1.5) 
vaccine_amount=[10,30,30,20,30,];%D:vaccination dosage
vaccine_time=[0,200,250,420,440];%T:vaccination time /days
vaccine_type=[2,2,2,2,2];%P: vaccination type(1:spike, 2:mRNA, 3:inactivated)

%----------------------

sub_varient=unique(vaccine_varient);
sub_varient=sub_varient(sub_varient~=0);
sub_X=X(sub_varient,sub_varient);
sub_self_neu=self_neu(sub_varient);


load('opt_parameter_human.mat')
parameter_set=parameter_set(1:22);%(model parameters）
[t,Rt,Agt,Abt,Ft,Mofft,Mont]=simIPD(vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,sub_X,sub_self_neu,parameter_set);%运行模拟
Mt=Mofft+Mont;

%----plot------
variant_color=[127,127,127;247,187,127;16,131,88;194,127,159;85,160,251;198,55,53]/255;
Yaxis_on=1;

f=figure(1);
set(f,'Position',[500,100,900,800])
position=[0.1,0.05,0.8,0.9];%position of the figure
dynaVac_pl(t,Agt,Abt,Ft,Mt,vaccine_varient,vaccine_time,variant_color,position,Yaxis_on)
