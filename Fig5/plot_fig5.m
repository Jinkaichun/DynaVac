clc,clear,close all
load('fig5_data.mat')
%%
close all
clc

%---------figure_position_parameters-------------
figure('color','w')
set(gcf,'Position',[50,50,1100,950])

w_left=0.41;
w_j=0.03;
w_dyna=(w_left-2*w_j)/3;

l_dyna=0.45;

l_b=0.05;
w_b=0.075;

l_scheme=0.05;
l_td_dyna=0.05;
l_b_up=0.04;
l_Ab_td=0.03;
l_td=(1-l_b_up-l_scheme-2*l_Ab_td-l_td_dyna-l_dyna-l_b)/3;


w_lr=0.1;
w_b_l=0.065;
w_right=1-w_b-w_b_l-w_lr-w_left;
l_ad=0.18;
l_inten=0.375;
l_i_run=0.07;
l_ad_run=0.075;
l_run=1-l_b-l_inten-l_i_run-l_ad_run-l_ad-l_b_up;
wj_run=0.03;
w_run_b=0;
w_right_run=1-w_b-w_run_b-w_lr-w_left;
w_run=(w_right_run-wj_run)/2;
w_add_in=0.005;
%------------------------

distance=(0:0.1:8);
ad_range=[0,8];
distancel=-2.^(-distance);
ad_rangel=-2.^(-ad_range);
boost_amount=(1:50);
variant_color=[127,127,127;204,51,17]/255;

lg2_titer_2=log2(titer_2_sample./titer_2_sample(:,1));
lg2_titer_1=log2(titer_1_sample./titer_1_sample(:,1));

[M_titer,I_titer]=min(lg2_titer_2,[],2);

num=0;
l_bound_titer_all=zeros(length(boost_amount),1);
r_bound_titer_all=zeros(length(boost_amount),1);
for check_amount=boost_amount
    num=num+1;
    I_titer_check=I_titer(boost_amount==check_amount);
    [~,l_bound_titer]=min(abs(lg2_titer_2(boost_amount==check_amount,1:I_titer_check)+1),[],2);
    [~,r_bound_titer]=min(abs(lg2_titer_2(boost_amount==check_amount,I_titer_check:end)+1),[],2);
    r_bound_titer=I_titer_check+r_bound_titer-1;
    l_bound_titer_all(num)=l_bound_titer;
    r_bound_titer_all(num)=r_bound_titer;
    
end

%-----------M-distance--------
check_amount=30;
ax_M=axes;
set(ax_M,'Position',[w_b,l_b+2*(l_td+l_Ab_td),w_left,l_td], 'Color', 'none')

plot(distancel,M_1_sample(boost_amount==check_amount,:),'LineWidth',1.5,'Color',variant_color(1,:))
hold on
plot(distancel,M_2_sample(boost_amount==check_amount,:),'LineWidth',1.5,'Color',variant_color(2,:))
ylabel('Memory B cell')

set(ax_M,'XScale','log')
xlim(ad_rangel)
ylim([0,max(M_1_sample(boost_amount==check_amount,:))])
% legend('primary','booster')
box off
Xlim =get(gca,'xlim');
Ylim=get(gca,'ylim');
line([Xlim(1) Xlim(2)],[Ylim(2) Ylim(2)],'color','k');
line([Xlim(2) Xlim(2)],[Ylim(1) Ylim(2)],'color','k');
ax_M.Color='none';
legend('PT specific','Variant specific','','')

%-----------Ab-distance--------
ax_Ab=axes;

set(ax_Ab,'Position',[w_b,l_b+1*(l_td+l_Ab_td),w_left,l_td], 'Color', 'none')

plot(distancel,Ab_1_sample(boost_amount==check_amount,:),'LineWidth',1.5,'Color',variant_color(1,:))
hold on
plot(distancel,Ab_2_sample(boost_amount==check_amount,:),'LineWidth',1.5,'Color',variant_color(2,:))
set(ax_Ab,'XScale','log')
ylabel('Antibody')
xlim(ad_rangel)

box off
Xlim =get(gca,'xlim');
Ylim=get(gca,'ylim');
line([Xlim(1) Xlim(2)],[Ylim(2) Ylim(2)],'color','k');
line([Xlim(2) Xlim(2)],[Ylim(1) Ylim(2)],'color','k');

legend('PT specific','Variant specific','','')
ax_Ab.Color='none';

%-----------titer-distance--------
ax_td=axes;
set(ax_td,'Position',[w_b,l_b,w_left,l_td], 'Color', 'none')

plot(distancel,lg2_titer_2(boost_amount==check_amount,:),'LineWidth',1.5,'Color',variant_color(2,:))
set(ax_td,'XScale','log')
hold on
%plot(distance,lg2_titer_1(pre_amount==check_amount,:),'LineWidth',1.5)
% scatter(distancel(I_titer(pre_amount==check_amount)),M_titer(pre_amount==check_amount),25,'r','filled')
% %plot([distance(I_titer),distance(I_titer)],[-0.65,M_titer],'r','LineStyle','--','LineWidth',1.5)
% text(distancel(I_titer(pre_amount==check_amount)),M_titer(pre_amount==check_amount)+0.3,['cross-neutralization=',num2str(100*round(2.^-distance(I_titer(pre_amount==check_amount)),2)),'%'],'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')

c_50pro=[ 0, 125, 150 ]/255;
plot(ad_rangel,[-1,-1],'LineStyle','--','LineWidth',1.5,'Color',c_50pro)
%text(ad_rangel(2),-1,'50% protection','Color',c_50pro,'VerticalAlignment','top')
xlim(ad_rangel)

ymin=-3;
ylim([ymin,0])


xlabel('Antigenic distance','FontSize',11,'FontWeight','normal')
ylabel('neutralization against variant')
box off
Xlim =get(gca,'xlim');
Ylim=get(gca,'ylim');
line([Xlim(1) Xlim(2)],[Ylim(2) Ylim(2)],'color','k');
line([Xlim(2) Xlim(2)],[Ylim(1) Ylim(2)],'color','k');

ax_td.Color='none';

l_bound_titer=l_bound_titer_all(boost_amount==check_amount);
r_bound_titer=r_bound_titer_all(boost_amount==check_amount);

w_prot=w_left*l_bound_titer/(ad_range(2)/0.1+1);
w_naive=w_left*(1-r_bound_titer/(ad_range(2)/0.1+1));
w_inval=w_left-w_prot-w_naive;

check_index=[11,68];
y_adj=0.2;
check_distance_all=[distance(check_index(1)),distance(I_titer_check),distance(check_index(2))];

% scatter(distancel(check_index(1)),lg2_titer_2(pre_amount==check_amount,check_index(1)),25,'filled','MarkerFaceColor',[0,0.5,0])
% 
% scatter(distancel(check_index(2)),lg2_titer_2(pre_amount==check_amount,check_index(2)),25,'b','filled')
% 
% text(distancel(check_index(1)),lg2_titer_2(pre_amount==check_amount,check_index(1))-y_adj,['cross-neutralization=',num2str(100*round(2.^-distance(check_index(1)),3)),'%'],'HorizontalAlignment','center',...
%     'VerticalAlignment','top')
% 
% text(distancel(check_index(2)),lg2_titer_2(pre_amount==check_amount,check_index(2))-y_adj,['cross-neutralization=',num2str(100*round(2.^-distance(check_index(2)),3)),'%'],'HorizontalAlignment','center',...
%     'VerticalAlignment','top')

y_tex=2300;

text(ax_Ab,-2^(-distance(l_bound_titer)/2),y_tex,{'Immune-Imprinting- ';'Protection Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.3,0],'FontSize',10,'FontWeight','bold')


text(ax_Ab,-2^(-(distance(l_bound_titer)+distance(r_bound_titer))/2),y_tex,{'Immune-Imprinting-';'Pitfall Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5,0,0],'FontSize',10,'FontWeight','bold')


text(ax_Ab,-2^(-(ad_range(2)+distance(r_bound_titer))/2),y_tex,{'Immune-Imprinting-';'Breakthrough Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.5],'FontSize',10,'FontWeight','bold')

% set(ax_td,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','50%','75%','90%','95%','97.5%','99%','99.5%','99.8%'})
% set(ax_Ab,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','50%','75%','90%','95%','97.5%','99%','99.5%','99.8%'})
% set(ax_M,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','50%','75%','90%','95%','97.5%','99%','99.5%','99.8%'})


set(ax_td,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','0.5','0.75','0.9','0.95','0.975','0.99','0.995','0.998'})
set(ax_Ab,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',[])
set(ax_M,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',[])

set(ax_td,'YTick',[-3,-2,-1,0],'YTickLabel',{'12.5%','25%','50%','100%'},'YTickMode','manual')
set(ax_td,'TickDir','out')
set(ax_Ab,'TickDir','out')
set(ax_M,'TickDir','out')

%-------background_protection--------
ax_pro=axes;
set(ax_pro,'Position',[w_b,l_b,w_prot,3*l_td+2*l_Ab_td])
pro_range=[0,distance(l_bound_titer)];

set(ax_pro, 'Color', 'none','XTick',[],'YTick',[],'XLim',pro_range,'Ylim',[0,1],'YColor','none');
ad_grad=(ad_range(1):0.1:distance(l_bound_titer));
y_grad=(0:0.1:1);
[AD,Y]=meshgrid(ad_grad,y_grad);
Z=ones(size(AD)).*lg2_titer_2(boost_amount==check_amount,1:l_bound_titer);
c_pro=[linspace(230,139,256)',linspace(225,189,256)',linspace(226,164,256)']/256;
colormap(ax_pro,c_pro);

ss=surface(AD,Y,Z,'EdgeColor','none','FaceAlpha',0.5,'FaceColor','interp');
uistack(ax_pro,'bottom')


%-------background_inval--------
ax_inval=axes;
set(ax_inval,'Position',[w_b+w_prot,l_b,w_inval,3*l_td+2*l_Ab_td])
inval_range=[distance(l_bound_titer+1),distance(r_bound_titer)];

set(ax_inval, 'Color', 'none','XTick',[],'YTick',[],'XLim',inval_range,'Ylim',[0,1]);
ad_grad=(distance(l_bound_titer+1):0.1:distance(r_bound_titer));
y_grad=(0:0.1:1);
[AD,Y]=meshgrid(ad_grad,y_grad);

Z=-ones(size(AD)).*lg2_titer_2(boost_amount==check_amount,l_bound_titer+1:r_bound_titer);

c_pit=[linspace(255,216,256)',linspace(230,121,256)',linspace(218,127,256)']/256;
colormap(ax_inval,c_pit);

surface(ax_inval,AD,Y,Z,'EdgeColor','none','FaceAlpha',0.4,'FaceColor','interp');
axis(ax_inval, 'off');
uistack(ax_inval,'bottom')



%-------background_naive--------
ax_naive=axes;
set(ax_naive,'Position',[w_b+w_prot+w_inval,l_b,w_naive,3*l_td+2*l_Ab_td])
naive_range=[distance(r_bound_titer+1),ad_range(2)];

set(ax_naive, 'Color', 'none','XTick',[],'YTick',[],'XLim',naive_range,'Ylim',[0,1]);
ad_grad=(distance(r_bound_titer+1):0.1:ad_range(2));
y_grad=(0:0.1:1);
[AD,Y]=meshgrid(ad_grad,y_grad);

Z=ones(size(AD)).*lg2_titer_2(boost_amount==check_amount,r_bound_titer+1:(ad_range(2)/0.1+1));

c_brk=[linspace(230,118,256)',linspace(256,157,256)',linspace(256,234,256)']/256;
colormap(ax_naive,c_brk);

surface(ax_naive,AD,Y,Z,'EdgeColor','none','FaceAlpha',0.35,'FaceColor','interp');
axis(ax_naive, 'off');
uistack(ax_naive,'bottom')

%-----------vaccine dynamics--------
vaccine_varient=[1,1,2,2];
vaccine_type=[2,2,2,2];
vaccine_amount=[3,3,30,30];
vaccine_time=[10,40,400,430]-10;

check_amount=30;

%check_distance_all=[1,4.2,7];

for num=(1:3)
    
    check_distance=check_distance_all(num);
    i=boost_amount==check_amount;
    j=distance==check_distance;

    t=t_amount_dis{i,j}-10;
    Agt=Ag_amount_dis{i,j};
    Abt=Ab_amount_dis{i,j};
    Ft=F_amount_dis{i,j};
    M_offt=M_off_amount_dis{i,j};
    M_ont=M_on_amount_dis{i,j};
    Mt=M_offt+M_ont;

    
    vaccine_correspond=[1;2];

    trunc_start=55;
    trunc_end=360;
    
    if num==1
        Yaxis_on=1;

    else
        Yaxis_on=0;

    end
    
    position=[w_b+(num-1)*(w_dyna+w_j),1-l_b_up-l_scheme-l_dyna,w_dyna,l_dyna];
    Ab_lim=[0,6700];
    dynaVac_pl_pos(t,Agt,Abt,Ft,Mt,vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,vaccine_correspond,variant_color,trunc_start,trunc_end,position,Yaxis_on,Ab_lim)
end

%-----------AD table--------

vaccines={'PT','Alpha/Beta','Delta','BA.1','BA.2/BA.5','XBB.1.5'};


pse={'PT','Alpha','Beta','Delta','BA.1','BA.2','BA.5','BF.7','BQ.1.1','XBB.1','XBB.1.5','JN.1'};
load('X_pse_yunda.mat')
%load('X_yunda.mat')
sub_pse=[(1:9),11,12];
X_pse=X_pse(sub_pse,:);
load('mix_parameters.mat')
X_pse(7:11,1)=X_pse(7:11,1).*parameter_set(25:29);
X_pse(10,2:4)=X_pse(10,2:4).*parameter_set(30:32)';

ad=-log2(X_pse);
ad(1,1)=0;
ad(2,2)=0;
ad(3,2)=0;
ad(4,3)=0;
ad(5,4)=0;
ad(6:7,5)=0;
ad(10,6)=0;

pse=pse(sub_pse);
ad=round(ad,2);

adl=1-2.^(-ad);

ax_ad=axes;
set(ax_ad,'Position',[1-w_b_l-w_right,1-l_b_up-l_ad,w_right,l_ad])
%axis off

%c_pro=othercolor('Greens3');
%c_pit=othercolor('Reds3');
%c_brk=othercolor('Bu_10');

ind_pro_ad=ad<distance(l_bound_titer+1);
ind_brk_ad=ad>distance(r_bound_titer);
%ind_pit=adl>=l_bound_titer_all|adl<=r_bound_titer_all;
min_pro=log2(0.5);
max_pro=0;
min_brk=log2(0.5);
max_brk=0;
max_pit=log2(0.5);
min_pit=min(lg2_titer_2(boost_amount==check_amount,:));
xx=(1:size(adl,2)+1)-1;
yy=flip((0:size(adl,1)));
[XX,YY]=meshgrid(xx,yy);
ZZ=ones(size(XX));
CC=ones([size(XX),3]);
for i=(1:size(adl,1))
    for j=(1:size(adl,2))
        dis_val=ad(i,j);
        [~,ind_val]=min(abs(dis_val-distance));
        val=lg2_titer_2(boost_amount==check_amount,ind_val);
        
        if ind_pro_ad(i,j)
            norm_value=(val-min_pro)/(max_pro-min_pro);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_pro(c_ind,:);
            
        elseif ind_brk_ad(i,j)
            norm_value=(val-min_brk)/(max_brk-min_brk);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_brk(c_ind,:);
        else
            norm_value=1-(val-min_pit)/(max_pit-min_pit);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_pit(c_ind,:);
        end
    end
end


surface(ax_ad,XX,YY,ZZ,CC,'FaceAlpha',0.35)

%h=heatmap(ad,'FontSize',0.00000001);

%colormap(h,(othercolor('OrRd8')));

%h.XData=vaccines;
%h.YData=pse;
%h.Colormap=cmp;
%h.Title='Antigenic distance';
%h.ColorbarVisible=0;

% ax_adtex=axes;
% set(ax_adtex,'Position',[1-w_b_l-w_right,1-l_b_up-l_ad,w_right,l_ad],'Color','none')
pse_fl=flip(pse);
xlim([0,6])
set(ax_ad,'XTick',(0.5:1:5.5),'XTickLabel',vaccines)
ylim([0,11])
set(ax_ad,'YTick',(0.5:1:10.5),'YTickLabel',pse_fl,'TickLength',[0;0])

[XX,YY]=meshgrid((0.5:1:5.5),(10.5:-1:0.5));

ad_str_array = arrayfun(@(x) num2str(round(x,4)), adl, 'UniformOutput', false);
text(ax_ad,XX(:),YY(:),ad_str_array(:),'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8)


pairs1=[2,3,4,5,6,10,11,11];
pairs2=[1,1,1,1,1,1,5,6];

ytex_ad=[-2.5,-2.5,-2.65,-2.5,-2.5,-2.45,-1.9,-2.2];
for i=(1:length(pairs1))
    ad_pair=ad(pairs1(i),pairs2(i));
    ad_pairl=-2.^(-ad_pair);
    scatter(ax_td, min(ad_pairl,ad_rangel(2)),ymin,'^','filled','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',0.2)   
    pair_tex=[vaccines{pairs2(i)},' vs. ',pse{pairs1(i)}];
    text(ax_td,min(ad_pairl,ad_rangel(2)),ytex_ad(i),pair_tex,'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0],'FontSize',7.5,'FontWeight','bold')
end

%-------Booster dosage--------

ax_inten=axes;
set(ax_inten,'Position',[1-w_b_l-w_right,w_b_l,w_right+w_add_in,l_inten])
[XX,YY]=meshgrid(distance,boost_amount);
XXl=-(2.^(-XX));
s=surface(XXl,YY,titer_2_sample,'EdgeColor','k','FaceColor','interp','EdgeAlpha',0.2);
%s=surface(XX,YY,lg2_titer_2,'EdgeColor','k','FaceColor','interp');
view(3)
set(ax_inten,'XScale','log')
set(ax_inten,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','0.5','0.75','0.9','0.95','0.975','0.99','0.995','0.998'})
set(ax_inten,'YDir','normal')
colormap(ax_inten,flip(othercolor('RdBu4')));


xlabel('Antigenic distance')
xlim(ad_rangel)
ylim([1,40])
ylabel('Booster dosage (\mug mRNA vaccine)')
zlabel('Combined neutralization against variant')

z_bound=-1200;
set(ax_inten,'ZTick',(0:500:3500),'ZTickLabel',(0:500:3500),'ZLim',[z_bound,3500])


CC=0.75*ones([size(XX),3]);


titer_2_sample_check=titer_2_sample(boost_amount==30,:);


ind_pro=distance<2.1;
ind_brk=distance>=5.5;
ind_pit=distance<5.5 &distance>=2.1;

min_pro=min(titer_2_sample_check(ind_pro));
max_pro=max(titer_2_sample_check(ind_pro));

min_pit=min(titer_2_sample_check(ind_pit));
max_pit=max(titer_2_sample_check(ind_pit));

min_brk=min(titer_2_sample_check(ind_brk));
max_brk=max(titer_2_sample_check(ind_brk));




for i=(1:length(distance))
    val=titer_2_sample_check(i);
    if ind_pro(i)
        norm_value=(val-min_pro)/(max_pro-min_pro);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC(:,i,:)=reshape(repmat(c_pro(c_ind,:),50,1),50,1,3);

        
    elseif ind_brk(i)
        norm_value=(val-min_brk)/(max_brk-min_brk);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC(:,i,:)=reshape(repmat(c_brk(c_ind,:),50,1),50,1,3);

        
    else
        norm_value=1-(val-min_pit)/(max_pit-min_pit);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC(:,i,:)=reshape(repmat(c_pit(c_ind,:),50,1),50,1,3);

    end
end


surface(ax_inten,XXl,YY,z_bound*ones(size(XX)),CC,'EdgeColor','none','FaceLighting','none','FaceAlpha',0.4)
hold on

c_50pro=[ 110, 110, 110 ]/255;

zz=(z_bound:100:3500);
[ZZ,YY2]=meshgrid(zz,boost_amount);

y_tex_inten=12.5;
z_tex=-900;
fsize=10;


%plot3(ax_inten,(-2^(-2.1))*ones(size(ZZ)),boost_amount,titer_2_sample(:,distance==2.1),'--','LineWidth',2,'Color',c_50pro);
%plot3(ax_inten,(-2^(-5.5))*ones(size(ZZ)),boost_amount,titer_2_sample(:,distance==5.5),'--','LineWidth',2,'Color',c_50pro);
surface(ax_inten,(-2^(-2.1))*ones(size(ZZ)),YY2,ZZ,'FaceColor',1.5*c_50pro,'FaceAlpha',0.5,'EdgeColor','none');
surface(ax_inten,(-2^(-5.5))*ones(size(ZZ)),YY2,ZZ,'FaceColor',1.5*c_50pro,'FaceAlpha',0.5,'EdgeColor','none');
text(ax_inten,-2^(-distance(22)/2),y_tex_inten,z_tex,{'Immune-Imprinting- ';'Protection Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.3,0],'FontSize',fsize,'FontWeight','bold')

text(ax_inten,-2^(-(distance(22)+distance(56))/2),y_tex_inten,z_tex,{'Immune-Imprinting-';'Pitfall Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5,0,0],'FontSize',fsize,'FontWeight','bold')

text(ax_inten,-2^(-(ad_range(2)+distance(56))/2),y_tex_inten,z_tex,{'Immune-Imprinting-';'Breakthrough Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.3],'FontSize',fsize,'FontWeight','bold')
plot3((-2^(-1))*ones(size(ZZ)),boost_amount,titer_2_sample(:,distance==1),'-','LineWidth',2,'Color',[0,0.5,0])
plot3((-2^(-4.3))*ones(size(ZZ)),boost_amount,titer_2_sample(:,distance==4.3),'-','LineWidth',2,'Color',[0.7,0,0])
plot3((-2^(-6.7))*ones(size(ZZ)),boost_amount,titer_2_sample(:,distance==6.7),'-','LineWidth',2,'Color',[0,0,0.5])

z_tex2=3500;
% text(ax_inten,-2^(-distance(22)/2),y_tex_inten,z_tex2,'Cross-Neutralization=50%','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.5,0],'FontSize',10,'FontWeight','bold')
% 
% text(ax_inten,-2^(-(distance(30)+distance(56))/2),y_tex_inten,z_tex2-2000,'Cross-Neutralization=5%','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.7,0,0],'FontSize',10,'FontWeight','bold')
% 
% text(ax_inten,-2^(-(ad_range(2)+distance(60))/2),y_tex_inten,z_tex2,'Cross-Neutralization=1%','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.5],'FontSize',10,'FontWeight','bold')

view(ax_inten,[-353.546197857206 13.59558772688]);

%-------MN_phase--------
ax_MN=axes;
set(ax_MN,'Position',[1-w_b_l-w_right,w_b_l+l_inten+l_i_run,w_run,l_run])
title(ax_MN,'1st boost (390 days)')
xlim([0,1])
ylim([0,0.4])

xlabel(ax_MN,'PT specific memory B cell level')
ylabel(ax_MN,'Variant specific naive B cell level')
set(ax_MN,'XTick',(0:0.2:1),'TickDir','out');
hold(ax_MN,'on')
uistack(ax_MN,'bottom')
box on
axis square
%-------MM_phase--------
ax_MM=axes;
set(ax_MM,'Position',[1-w_run_b-w_run,w_b_l+l_inten+l_i_run,w_run,l_run])


hold(ax_MM,'on')
box on

immune_span=9;






check_index=[16,51,61];
check_distance=distance(check_index);
distance_c=distance(15:61);
titer_2_sample_check=titer_2_sample_check(15:61);

ind_pro=distance_c<2.1 ;
ind_brk=distance_c>=5.5; 
ind_pit=distance_c<5.5 &distance_c>=2.1;

min_pro=min(titer_2_sample_check(ind_pro));
max_pro=max(titer_2_sample_check(ind_pro));

min_pit=min(titer_2_sample_check(ind_pit));
max_pit=max(titer_2_sample_check(ind_pit));

min_brk=min(titer_2_sample_check(ind_brk));
max_brk=max(titer_2_sample_check(ind_brk));


CC_curve=zeros(length(distance_c),3);


for i=(1:length(distance_c))
    val=titer_2_sample_check(i);
    if ind_pro(i)
        norm_value=(val-min_pro)/(max_pro-min_pro);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC_curve(i,:)=c_pro(c_ind,:);
        
    elseif ind_brk(i)
        norm_value=(val-min_brk)/(max_brk-min_brk);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end

        CC_curve(i,:)=c_brk(c_ind,:);
        
    else
        norm_value=1-(val-min_pit)/(max_pit-min_pit);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC_curve(i,:)=c_pit(c_ind,:);
    end
end

CC_curve=0.95*CC_curve;
i=boost_amount==30;


%check_index(2)=61;
%check_distance=[distance(check_index(1)),distance(I_titer_check),distance(check_index(2))];

% CC_curve=CC_curve(check_index(1):check_index(end),:);

cmp=colormap(ax_MN,CC_curve);
colormap(ax_MM,CC_curve);
dis_cmp=cmp;
%dis_cmp=cmp(round(linspace(1,length(CC_curve),check_index(end)-check_index(1)+1)),:);


num=0;
num_plot=0;
dis_max_M_bos_1=[];
dis_max_F_bos_1=[];
dis_s_pre_M_bos_2=[];
dis_s_new_M_bos_2=[];
dis_f_pre_M_bos_2=[];
dis_f_new_M_bos_2=[];

for dis=distance_c
    num=num+1;
    j=distance==dis;
    
    t=t_amount_dis{i,j};
    Agt=Ag_amount_dis{i,j};
    Abt=Ab_amount_dis{i,j};
    Ft=F_amount_dis{i,j};
    M_offt=M_off_amount_dis{i,j};
    M_ont=M_on_amount_dis{i,j};
    Mt=M_offt+M_ont;

    t_bos_1=t(t>400&t<400+immune_span); %第一次boost免疫反应时间段
    t_bos_2=t(t>430&t<430+immune_span); %第二次boost免疫反应时间段
   
    Mt_bos_1=Mt(:,t>400&t<400+immune_span);
    Mt_bos_2=Mt(:,t>430&t<430+immune_span);
    
    pre_Mt_bos_1=Mt_bos_1(1,:);
    pre_Mt_bos_2=Mt_bos_2(1,:);
    new_Mt_bos_2=Mt_bos_2(2,:);
    
    [~,ind]=max(pre_Mt_bos_2);
    pre_Mt_bos_2=pre_Mt_bos_2(1:ind);
    new_Mt_bos_2=new_Mt_bos_2(1:ind);
    
    t_bos_2=t_bos_2(1:ind);
    
    
    Ft_bos_1=Ft(:,t>400&t<400+immune_span);
    new_Ft_bos_1=Ft_bos_1(2,:);
    
    [max_F_bos_1,max_F_id]=max(new_Ft_bos_1);
    max_M_bos_1=pre_Mt_bos_1(max_F_id);
    
    dis_s_pre_M_bos_2=[dis_s_pre_M_bos_2,pre_Mt_bos_2(1)];
    dis_s_new_M_bos_2=[dis_s_new_M_bos_2,new_Mt_bos_2(1)];
    
    dis_f_pre_M_bos_2=[dis_f_pre_M_bos_2,pre_Mt_bos_2(end)];
    dis_f_new_M_bos_2=[dis_f_new_M_bos_2,new_Mt_bos_2(end)];
    
    dis_max_M_bos_1=[dis_max_M_bos_1,max_M_bos_1];
    dis_max_F_bos_1=[dis_max_F_bos_1,max_F_bos_1];   
    
    if sum(dis==check_distance)
          
%         scale=0.025;
%         quiver_width=1.1;
%         int=150;
        start=140;
        num_plot=num_plot+1;

        plot(ax_MN,pre_Mt_bos_1,new_Ft_bos_1,'Color',dis_cmp(num,:))
        d_pre_Mt_bos_1=diff(pre_Mt_bos_1);
        d_new_Ft_bos_1=diff(new_Ft_bos_1);
        scatter(ax_MN,max_M_bos_1,max_F_bos_1,20,'filled','MarkerFaceColor',dis_cmp(num,:))
        

        scale=0.025;
        quiver_width=1;
        int=100;

        plot(ax_MM,pre_Mt_bos_2,new_Mt_bos_2,'Color',dis_cmp(num,:))

        d_pre_Mt_bos_2=diff(pre_Mt_bos_2);
        d_new_Mt_bos_2=diff(new_Mt_bos_2);
       
        scatter(ax_MM,pre_Mt_bos_2(1),new_Mt_bos_2(1),20,'x','MarkerFaceColor',dis_cmp(num,:),'MarkerEdgeColor',dis_cmp(num,:))
        scatter(ax_MM,pre_Mt_bos_2(end),new_Mt_bos_2(end),20,'s','MarkerFaceColor',dis_cmp(num,:),'MarkerEdgeColor',dis_cmp(num,:))
        t_bos_1=t_bos_1-10;
        t_bos_2=t_bos_2-10;
%         figure(2)
%         
%         set(gcf,'Position',[50,50,600,950])
%         subplot(3,2,2*num_plot-1)
%         yyaxis left
%         plot(t_bos_1,pre_Mt_bos_1,'Color',variant_color(1,:),'LineWidth',1.5)
%         ylim([0,1])
%         ylabel('PT specific memory B cell level')
%         set(gca,'YColor',variant_color(1,:))
%         yyaxis right
%         plot(t_bos_1,new_Ft_bos_1,'Color',variant_color(2,:),'LineWidth',1.5)
%         ylim([0,0.4])
%         ylabel('Variant specific naive B cell level')
%         xlabel('Time (Days)')
%         set(gca,'YColor',variant_color(2,:))
%         xlim([min(t_bos_1),max(t_bos_1)])
%         set(gca,'TickDir','out')
%         subplot(3,2,2*num_plot)
%         yyaxis left
%         plot(t_bos_2,pre_Mt_bos_2,'Color',variant_color(1,:),'LineWidth',1.5)
%         ylim([0,1])
%         ylabel('PT specific memory B cell level')
%         set(gca,'YColor',variant_color(1,:))
%         yyaxis right
%         plot(t_bos_2,new_Mt_bos_2,'Color',variant_color(2,:),'LineWidth',1.5)
%         ylim([0,0.4])
%         ylabel('Variant specific memory B cell level')
%         xlabel('Time (Days)')
%         set(gca,'YColor',variant_color(2,:))
%         xlim([min(t_bos_2),max(t_bos_2)])
%         set(gca,'TickDir','out')
    end
end

dis_max_F_bos_1(end)=NaN;

c=(1:length(dis_max_F_bos_1));
patch(ax_MN,dis_max_M_bos_1,dis_max_F_bos_1,c,'EdgeColor','interp','Marker','none','MarkerFaceColor','flat','MarkerSize',0.01,'LineStyle',':','LineWidth',1);
l1=colorbar(ax_MN);
l1.Ticks=[];
%l1.TicksLabel=[5,5,5,5,5];
%l1.Label.String='Antigenic distance';

title(ax_MM,'2nd boost (420 days)')
xlim(ax_MM,[0,1])
ylim(ax_MM,[0,0.4])

xlabel(ax_MM,'PT specific memory B cell level')
ylabel(ax_MM,'Variant specific memory B cell level')
set(ax_MM,'XTick',(0:0.2:1),'TickDir','out');

dis_s_new_M_bos_2(end)=NaN;
c=(1:length(dis_s_new_M_bos_2));
patch(ax_MM,dis_s_pre_M_bos_2,dis_s_new_M_bos_2,c,'EdgeColor','interp','Marker','none','MarkerFaceColor','flat','MarkerSize',0.01,'LineStyle',':','LineWidth',1);

dis_f_new_M_bos_2(end)=NaN;
c=(1:length(dis_f_new_M_bos_2));

patch(ax_MM,dis_f_pre_M_bos_2,dis_f_new_M_bos_2,c,'EdgeColor','interp','Marker','none','MarkerFaceColor','flat','MarkerSize',0.01,'LineStyle','-.','LineWidth',1);


l2=colorbar(ax_MM);
l2.Ticks=[];
l2.Label.String='Antigenic distance';
%axis square

%---------text------------
letter_size=10;
ax_tex=axes;
set(ax_tex,'Position',[0,0,1,1])
axis off
% text(w_b/3,1-l_b_up/2,'A','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% text(w_b/3,1-l_b_up-l_td-l_Ab_td/2,'B','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% text(w_b/3,1-l_b_up-2*l_td-l_Ab_td-l_Ab_td/2,'C','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% 
% text(w_b/2,1-l_b_up-3*l_td-2*l_Ab_td-l_td_dyna,'D','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% text(w_b+w_dyna+w_j/2,1-l_b_up-3*l_td-2*l_Ab_td-l_td_dyna,'E','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% text(w_b+2*w_dyna+3*w_j/2,1-l_b_up-3*l_td-2*l_Ab_td-l_td_dyna,'F','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% 
% 
% text(w_b+w_left+w_lr/1.5,1-l_b_up/2,'G','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% 
% 
% text(w_b+w_left+w_lr/1.5,1-l_b_up-l_ad-l_ad_run/2,'H','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% 
% text(w_b+w_left+w_lr+w_run+wj_run/2,1-l_b_up-l_ad-l_ad_run/2,'I','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
% 
% text(w_b+w_left+w_lr/1.5,1-l_b_up-l_ad-l_ad_run-l_run-l_i_run/2,'J','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')

%text(w_b+w_left/2,l_b/2,'Time (days)','FontSize',12,'FontWeight','normal','HorizontalAlignment','center','VerticalAlignment','middle')

%text(w_b+w_left+w_lr+w_right/2,1-l_b_up/2,'Antigenic distance','FontSize',10,'Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle')