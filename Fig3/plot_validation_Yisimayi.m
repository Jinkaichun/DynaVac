clc,clear,close all

%-----load parameters----- 
load('cyl_parameters.mat')

%-----input batch vaccination-------
load('batch_information_all.mat')
load('titers_all.mat')

load('X.mat')
self_neu=[1,1,0.654107249,0.253360062,1]';

parameter_set(21)=1;
%parameter_set(33)=0.4*parameter_set(20);
parameter_set_new=parameter_set(1:22);

self_neu(5)=parameter_set(33);
%parameter_set_new(19)=0.05;

% load('ga_opt_parameters_ode_broad.mat')
% parameter_set_new=parameter_reduction(parameter_set);

vaccine_correspond=(1:5)';
X_pse=X(1:4,:);

%----run-----
[batch_t,batch_Rt,batch_Agt,batch_Abt,batch_Ft,batch_Mofft,batch_Mont]=batch_simIPD(batch_vaccine_varient,batch_vaccine_time,batch_vaccine_type,batch_vaccine_amount,X,self_neu,parameter_set_new,vaccine_correspond);
[batch_titer_lgFC,batch_titer_est_lgFC,titer_group,titer_sample]=batch_getFC_cyl(batch_vaccine_varient,batch_vaccine_time,batch_t,batch_Abt,batch_sample,titers,X_pse,self_neu,vaccine_correspond);



%% plots

n_pse=size(X_pse,1);


l_b_up=0.075;
l_b=0.02;


w_b_r=0.1;

w_b=0.1;
l_jj_pse=0.055;
l_comb=0.2;
l_leg=0.2;
l_jj_leg=0.01;


l_pse=(1-l_comb-l_b_up-l_b-(n_pse-1)*l_jj_pse-l_jj_pse-l_leg-l_jj_leg)/n_pse;

l_pse_all=1-l_comb-l_b_up-l_b-l_jj_pse-l_leg-l_jj_leg;

w_left=1-w_b_r-w_b;


figure(1)
set(gcf,'Position',[50,50,600,850])
set(gcf,'Color','w')
n_titer=size(batch_titer_lgFC,2);

order_index=order_vaccination(batch_vaccine_varient,batch_vaccine_time,batch_vaccine_type,batch_vaccine_amount,titer_group,titer_sample);

labels={'D614G','BA.5','BQ.1.1','XBB.1.5'};


c_exp=[255 165 0]/255;
% c_est=[107 62 151]/255;
c_est=[0 0 139]/255;
alpha=0.3;
msize=50;
for i=(1:4)
    
    ax_pse(i)=axes;
    set(ax_pse(i),'Position',[w_b,1-l_b_up-l_pse-(i-1)*(l_pse+l_jj_pse),w_left,l_pse]) 
    hold on
    scatter((1:n_titer),batch_titer_lgFC(i,order_index),msize,'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_exp,'MarkerEdgeColor',c_exp,'LineWidth',1)
    scatter((1:n_titer),batch_titer_est_lgFC(i,order_index),msize,'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_est,'MarkerEdgeColor',c_est,'LineWidth',1)
    plot((1:n_titer),batch_titer_lgFC(i,order_index),'-','Color',c_exp)
    plot((1:n_titer),batch_titer_est_lgFC(i,order_index),'-','Color',c_est)
    xlabel(labels(i))
    ylim([-9,8])
    xlim([0.5,n_titer+0.5])
    set(gca,'XTick',(1:n_titer))
   
end
lgd=legend('Experiment','Model');
set(lgd,'Parent',gcf)
set(lgd,'Units','normalized')
set(lgd,'Position',[0.78,0.945,0.04 0.04])
%-------vac_combo----------
ax_comb=axes;
set(ax_comb,'Position',[w_b,l_b+l_leg+l_jj_leg,w_left,l_comb])
hold on
varient_color=[0.5,0.5,0.5;1,0,0;0,1,0;0,0,1;0.5,0.5,0.5];

num=0;
scatter_sizes=[12,35,75]; %对应1ug 3ug 10ug
for i=order_index
    num=num+1;
    group=titer_group(i);
    sample=titer_sample(i);
    
    vaccine_varient_i= batch_vaccine_varient(group,1:sample+1);
    vaccine_time_i= batch_vaccine_time(group,1:sample+1);
    vaccine_type_i= batch_vaccine_type(group,1:sample+1);
    vaccine_amount_i= batch_vaccine_amount(group,1:sample+1);
    vaccine_interval=vaccine_time_i(2:end)-vaccine_time_i(1:end-1);
    interval_str=cell(1,sample);
    size_i=scatter_sizes(1).*(vaccine_amount_i==1)+scatter_sizes(2).*(vaccine_amount_i==3)+scatter_sizes(3).*(vaccine_amount_i==10);
    
    type_pr=(vaccine_type_i==1);
    type_mrna=(vaccine_type_i==2);
    
    for j=(1:sample)
        interval_str{j}=num2str(vaccine_interval(j));
    end
    XX=num*ones(1,sample+1);
    YY=(70:-20:70-20*sample);
    SZ=size_i;
    CC=varient_color(vaccine_varient_i,:);
    scatter(XX(type_pr),YY(type_pr),SZ(type_pr),CC(type_pr,:),'filled')
    hold on
    scatter(XX(type_mrna),YY(type_mrna),SZ(type_mrna),CC(type_mrna,:),'s','filled')

    text(num*ones(1,sample),(60:-20:80-20*sample),interval_str,'FontSize',6.5,'HorizontalAlignment','center')

end


xlim([0.5,n_titer+0.5])
ylim([0,80]);
set(gca,'YTick',[])
set(gca,'XTick',(1:n_titer))
ylabel('vaccine combo')

%---------legend----------
ax_leg=axes;
set(ax_leg,'Position',[w_b,l_b,w_left,l_leg])

vaccine_inform{1}='CoronaVac';
vaccine_inform{2}='BA.5 spike';
vaccine_inform{3}='BQ.1.1 spike';
vaccine_inform{4}='XBB.1.5 spike';

vaccine_inform{5}='WT spike mRNA';
vaccine_inform{6}='BA.5 spike mRNA';
vaccine_inform{7}='BQ.1.1 spike mRNA';
vaccine_inform{8}='XBB.1.5 spike mRNA';

amount{1}='1 ug';
amount{2}='3 ug';
amount{3}='10 ug';

ybias=0.4;
scatter(10*ones(1,4),linspace(20,5,4),scatter_sizes(2),varient_color(1:4,:),'filled')
hold on
scatter(38*ones(1,4),linspace(20,5,4),scatter_sizes(2),varient_color(1:4,:),'s','filled')
%scatter(10*ones(1,4),(35:-10:5),scatter_sizes(2),varient_color,'filled')
text(10*ones(1,4)+3,linspace(20,5,4)+ybias,vaccine_inform(1:4),'FontWeight','bold')
text(38*ones(1,4)+3,linspace(20,5,4)+ybias,vaccine_inform(5:8),'FontWeight','bold')

scatter(77*ones(1,3),linspace(20,5,3),scatter_sizes,varient_color(1,:),'filled')
scatter(77*ones(1,3)+3,linspace(20,5,3),scatter_sizes,varient_color(1,:),'s','filled')
text(77*ones(1,3)+6,linspace(20,5,3)+ybias,amount,'FontWeight','bold')

axis equal
xlim([5,95])
ylim([5,20])

axis off

ax_text=axes;
set(ax_text,'Position',[0,0,1,1])
axis off

text(w_b/2,1-l_pse_all/2-l_b_up,'log_2(titer FC)','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12,'FontWeight','normal','Rotation',90)

figure(2)
batch_dataset=zeros(size(batch_titer_lgFC));
draw_R_square(gca,batch_titer_lgFC,batch_titer_est_lgFC,batch_dataset)
%save('lgFC_yunda_pa_cyl.mat','batch_titer_lgFC','batch_titer_est_lgFC','batch_dataset')