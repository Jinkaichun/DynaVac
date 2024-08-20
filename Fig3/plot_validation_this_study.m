clc,clear,close all

%-----load parameters----- 
load('yd_parameters.mat')

parameter_set(21)=1;

parameter_set_new=parameter_set(1:22);
parameter_set_new([2,20])=parameter_set(23:24);

%-----input batch vaccination-------
load('batch_information_yunda.mat')
load('X_yunda.mat')
load('self_neu_yunda.mat');
self_neu=self_neu';
load('vaccine_correspond_yunda.mat')
load('titers_yunda.mat')
load('X_pse_yunda.mat')


pse={'WT','Alpha','Beta','Delta','BA.1','BA.2','BA.5','BF.7','BQ.1.1','XBB.1','XBB.1.5','JN.1'};
sub_pse=[(1:9),11,12];
X_pse=X_pse(sub_pse,:);
pse=pse(sub_pse);

X(5,1)=X(5,1)*sqrt(parameter_set(25));
X(6,1)=X(6,1)*parameter_set(28);
X(6,2)=X(6,2)*parameter_set(30);
X(6,3)=X(6,3)*parameter_set(31);
X(6,4)=X(6,4)*parameter_set(32);

X_pse(7:11,1)=X_pse(7:11,1).*parameter_set(25:29);
X_pse(10,2:4)=X_pse(10,2:4).*parameter_set(30:32)';

titers{1}([7:9,11,12])=titers{1}([7:9,11,12]).*parameter_set(25:29);
titers{2}(11)=titers{2}(11)*parameter_set(30);
titers{3}(11)=titers{3}(11)*parameter_set(31);
titers{4}(11)=titers{4}(11)*parameter_set(32);

%----run-----
[batch_t,batch_Rt,batch_Agt,batch_Abt,batch_Ft,batch_Mofft,batch_Mont]=batch_simIPD(batch_vaccine_varient,batch_vaccine_time,batch_vaccine_type,batch_vaccine_amount,X,self_neu,parameter_set_new,vaccine_correspond);
[batch_titer_lgFC,batch_titer_est_lgFC,titer_group,titer_time,titer_sample]=batch_getFC_yd(batch_vaccine_varient,batch_vaccine_time,batch_t,batch_Abt,batch_sample,titers,X_pse,self_neu,vaccine_correspond,sub_pse);


%%
clc,close all
figure(1)
set(gcf,'Position',[50,50,600,850])
set(gcf,'Color','w')
f_pos=get(gcf,'Position');
n_pse=size(X_pse,1);


l_b_up=0.03;
l_b=0.075;
w_b_r=0.175;


w_b=0.1;
l_jj_pse=0.03;
l_comb=0.2;

l_pse=(1-l_comb-l_b_up-l_b-(n_pse-1)*l_jj_pse-l_jj_pse)/n_pse;
w_leg=0.05;
w_jj_leg=0.035;
l_b_leg=0.075;
w_left=1-w_b_r-w_b;

w_jj_leg_pse=0.01;
l_leg_pse=0.05;


%----------%

n_titer=size(batch_titer_lgFC,2);
labels=pse;

order_index=[1,2,9,10,3,11,12,4,13,14,5,15,16,8,21,22,6,17,18,7,19,20,23,25,26,24,27,28];
order_index_3025=order_index(17:19);
order_index_3033=order_index(14:16);
% 
% order_index(14:16)=order_index_3025;
% order_index(17:19)=order_index_3033;


alpha=0.3;

c_exp=[255 165 0]/255;
% c_est=[107 62 151]/255;
c_est=[0 0 139]/255;

for i=(1:size(X_pse,1))
    ax_pse(i)=axes;
    set(ax_pse(i),'Position',[w_b,1-l_b_up-l_pse-(i-1)*(l_pse+l_jj_pse),w_left,l_pse]) 
    hold on
    
    scatter((1:n_titer),batch_titer_lgFC(i,order_index),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_exp,'MarkerEdgeColor',c_exp,'LineWidth',1)
    scatter((1:n_titer),batch_titer_est_lgFC(i,order_index),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_est,'MarkerEdgeColor',c_est,'LineWidth',1)
    plot((1:n_titer),batch_titer_lgFC(i,order_index),'-','Color',c_exp)
    plot((1:n_titer),batch_titer_est_lgFC(i,order_index),'-','Color',c_est)
    
    
    xlabel(labels(i),'FontSize',10)
     
    
    min_FC=min([batch_titer_lgFC(i,:),batch_titer_est_lgFC(i,:)],[],2,'includenan');
    max_FC=max([batch_titer_lgFC(i,:),batch_titer_est_lgFC(i,:)],[],2,'includenan');

    min_FC=min(min_FC);
    max_FC=max(max_FC);
    
    min_FC=min(min_FC,-5);
    max_FC=max(max_FC,5);
    
    lgFC_width=[min_FC-1,max_FC+1];
    ylim(lgFC_width)

    xlim([0.5,n_titer+0.5])
    set(gca,'XTick',(1:n_titer))
    ax=gca;
    ax.FontSize=5;
    ax.XLabel.FontSize=8;

    if i==6
        ylabel('log_2(titer FC)','FontSize',12)
    end
end
% lgd=legend('experiment','prediction');
% set(lgd,'Parent',gcf)
% set(lgd,'Units','normalized')
% set(lgd,'Position',[0.78,0.945,0.04 0.04])

ax_comb=axes;
set(ax_comb,'Position',[w_b,l_b,w_left,l_comb])
hold on
varient_color=[127,127,127;247,187,127;16,131,88;194,127,159;85,160,251;198,55,53]/255;


num=0;
y_space=4;

for i=order_index
    num=num+1;
    group=titer_group(i);
    sample=titer_sample(i);
    vaccine_varient_i= batch_vaccine_varient(group,1:sample);
    vaccine_varinet_i=vaccine_correspond(vaccine_varient_i,:)';
    
    
    vaccine_time_i= batch_vaccine_time(group,1:sample);
    vaccine_type_i= batch_vaccine_type(group,1:sample);
    vaccine_amount_i= batch_vaccine_amount(group,1:sample);
    vaccine_interval=vaccine_time_i(2:end)-vaccine_time_i(1:end-1);
    interval_str=cell(1,sample-1);

    
    for j=(1:sample-1)
        interval_str{j}=num2str(vaccine_interval(j));
    end
    
    for k=(1:sample)
        injection=vaccine_varinet_i(:,k);
        injection=injection(injection~=0);
        XX=num*ones(size(injection));
        valent=length(injection);
        y_center=110-20*k;
        YY=(y_center-(valent-1)*y_space/2:y_space:y_center+(valent-1)*y_space/2);
        CC=varient_color(injection,:);
        if injection==1
            scatter(XX,YY,50,CC,'filled')
        else
            scatter(XX,YY,80,CC,'s','filled')
        end
        if valent==2
            y_white=5.5;
            YY=y_center+[-y_white,y_white];
            scatter(XX,YY,80,'s','filled','MarkerFaceColor','w','MarkerEdgeColor','none')
        end
        
        
    end
        
    text(num*ones(1,sample-1),(80:-20:120-20*sample),interval_str,'FontSize',8,'HorizontalAlignment','center')
end

xlim([0.5,n_titer+0.5])
ylim([0,100]);
set(gca,'YTick',[])
set(gca,'XTick',(1:n_titer))
xlabel('vaccine combo')


%--------legend_pse--------
ax_legend=axes;
set(ax_legend,'Position',[w_b+w_left+w_jj_leg_pse,1-l_b_up-l_leg_pse,w_leg,l_leg_pse])
hold on
scatter(5,8,'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_exp,'MarkerEdgeColor',c_exp,'LineWidth',1)
scatter(5,5,'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_est,'MarkerEdgeColor',c_est,'LineWidth',1)
axis off
text(6.5,8,'Experiment','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
text(6.5,5,'Model','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
xlim([4,10])
ylim([4,11])

%--------legend_vaccine--------
ax_legend=axes;
set(ax_legend,'Position',[w_b+w_left+w_jj_leg,l_b_leg,w_leg,l_comb])
vaccine_inform{1}='CoronaVac';
vaccine_inform{2}='Alpha/Beta mRNA';
vaccine_inform{3}='Delta mRNA';
vaccine_inform{4}='BA.1 mRNA';
vaccine_inform{5}='BA.2/4/5 mRNA';


vaccine_inform{1}='CoronaVac';
vaccine_inform{2}='RQ3013';
vaccine_inform{3}='RQ3014';
vaccine_inform{4}='RQ3021';
vaccine_inform{5}='RQ3019';
vaccine_inform{6}='RQ3033';

vaccine_inform{7}='RQ3025';
vaccine_inform{8}='RQ3027';

ybias=0.4;

y_space=0.61;
y_white=1.4;
y=linspace(30,5,8);
hold on
scatter(10*ones(1,1),y(1),50,varient_color(1,:),'filled')
scatter(10*ones(1,5),y(2:6),80,varient_color(2:6,:),'s','filled')

scatter(10*ones(1,2),y(7:8)+y_space,80,varient_color(5:6,:),'s','filled')
scatter(10*ones(1,2),y(7:8)-y_space,80,varient_color([2,2],:),'s','filled')


scatter(10*ones(1,2),y(7:8)+y_white,80,'s','filled','MarkerFaceColor','w','MarkerEdgeColor','none')
scatter(10*ones(1,2),y(7:8)-y_white,80,'s','filled','MarkerFaceColor','w','MarkerEdgeColor','none')

text(10*ones(1,8)+3,linspace(30,5,8)+ybias,vaccine_inform(1:8),'FontWeight','normal','FontSize',8)

%axis equal
xlim([10,20])
ylim([3,31])

axis off
figure(2)

titers{1}([7:9,11,12])=NaN;
titers{2}(11)=NaN;
titers{3}(11)=NaN;
titers{4}(11)=NaN;

[batch_titer_lgFC,batch_titer_est_lgFC]=batch_getFC_yd(batch_vaccine_varient,batch_vaccine_time,batch_t,batch_Abt,batch_sample,titers,X_pse,self_neu,vaccine_correspond,sub_pse);

batch_dataset=ones(size(batch_titer_lgFC));
draw_R_square(gca,batch_titer_lgFC,batch_titer_est_lgFC,batch_dataset)
%save('lgFC_cyl_pa_yunda.mat','batch_titer_lgFC','batch_titer_est_lgFC','batch_dataset')