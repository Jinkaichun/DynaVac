clc,clear,close all
load('fig6b_data.mat')
%%
interval=(1:36)*30;
distance=(0:0.1:8);

interval_range=[interval(1),interval(end)];

lg2_titer_2=log2(titer_2_sample./titer_2_sample(:,1));
lg2_titer_1=log2(titer_1_sample./titer_1_sample(:,1));

[M_titer,I_titer]=min(lg2_titer_2,[],2);

num=0;
l_bound_titer_all=zeros(length(interval),1);
r_bound_titer_all=zeros(length(interval),1);

for check_interval=interval
    num=num+1;
    I_titer_check=I_titer(interval==check_interval);
    [~,l_bound_titer]=min(abs(lg2_titer_2(interval==check_interval,1:I_titer_check)+1),[],2);
    [~,r_bound_titer]=min(abs(lg2_titer_2(interval==check_interval,I_titer_check:end)+1),[],2);
    r_bound_titer=I_titer_check+r_bound_titer-1;
    l_bound_titer_all(num)=l_bound_titer;
    r_bound_titer_all(num)=r_bound_titer;
    
end

clc,close all

distancel=-2.^(-distance);
ad_range=[0,8];
ad_rangel=-2.^(-ad_range);
figure(2)
set(gcf,'Position',[400,200,900,750])
ax_inten=axes;
set(ax_inten,'Position',[0.15,0.1,0.7,0.75])

[XX,YY]=meshgrid(distance,interval);
XXl=-(2.^(-XX));
s=surface(XXl,YY,lg2_titer_2,'EdgeColor','k','FaceColor','interp','EdgeAlpha',0.2);
%s=surface(XX,YY,lg2_titer_2,'EdgeColor','k','FaceColor','interp');
view(3)
set(ax_inten,'XScale','log')
set(ax_inten,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','0.5','0.75','0.9','0.95','0.975','0.99','0.995','0.998'})
set(ax_inten,'YDir','normal')
colormap(ax_inten,flip(othercolor('RdBu4')));
% s.CData=map;

xlabel('Antigenic distance')
xlim(ad_rangel)
ylim(interval_range)
ylabel('Interval between primary and booster vaccination \it v \rm (months)')
zlabel('Combined neutralization to variant (relative)')

c_50pro=[ 150, 150, 150 ]/255;
% [XXX,YYY]=meshgrid(ad_range,interval_range);
% XXXl=-(2.^(-XXX));


c_interval_3=1.3*[0.4940 0.1840 0.5560];
c_interval_2=0.9*[0.4940 0.1840 0.5560];
c_interval_1=0.75*[0.4940 0.1840 0.5560];

hold on
plot3(distancel,90*ones(size(distancel)),lg2_titer_2(interval==90,:),'-','LineWidth',2,'Color',c_interval_1)
plot3(distancel,360*ones(size(distancel)),lg2_titer_2(interval==360,:),'-','LineWidth',2,'Color',c_interval_2)
plot3(distancel,720*ones(size(distancel)),lg2_titer_2(interval==720,:),'-','LineWidth',2,'Color',c_interval_3)




text(ax_inten,-1,90,0,'\it v=3','VerticalAlignment','bottom','HorizontalAlignment','right','Color',c_interval_1,'FontWeight','bold','FontSize',10)
text(ax_inten,-1,360,0,'\it v=12','VerticalAlignment','bottom','HorizontalAlignment','right','Color',c_interval_2,'FontWeight','bold','FontSize',10)
text(ax_inten,-1,720,0,'\it v=24','VerticalAlignment','bottom','HorizontalAlignment','right','Color',c_interval_3,'FontWeight','bold','FontSize',10)





%text(ax_inten,0,30,'Booster interval=3 months','VerticalAlignment','bottom','HorizontalAlignment','left','Color',c_amount,'FontWeight','bold')

z_bound=-8;
set(ax_inten,'ZTick',[log2(0.025),log2(0.05),log2(0.1),-2,-1,0],'ZTickLabel',{'2.5%','5%','10%','25%','50%','100%'},'ZLim',[z_bound,1])
set(gca,'YTick',(90:90:interval_range(2)),'YTickLabel',(3:3:36),'YLim',[2,26]*30)

x_q=(0:0.05:8);
y_q=(1:0.1:36)*30;
[XX,YY]=meshgrid(x_q,y_q);
XXl=-(2.^(-XX));

%-----plot left bound------
xx_left=distance(l_bound_titer_all);
yy_left=interval;
xx_left_q = interp1(yy_left,xx_left,y_q);

xx_left_ql=-2.^(-xx_left_q);

% xx_left_all=distance(l_bound_titer_all);
% 
% xx_left=distance(l_bound_titer_all);
% yy = smooth(xx_left,interval,'sgolay');
% xx_leftl=-2.^(-xx_left);

% xx_left=[xx_left,xx_left(end)];

%yy=[yy;40];

plot3(xx_left_ql,y_q,z_bound*ones(size(y_q)),'--','LineWidth',2,'Color',c_50pro)

%-----plot right bound------


xx_right=distance(r_bound_titer_all);
yy_right=interval;
xx_right_q = interp1(yy_right,xx_right,y_q);

xx_right_ql=-2.^(-xx_right_q);

xx_right=distance(r_bound_titer_all);
xx_rightl=-2.^(-xx_right);
yy = smooth(xx_rightl,interval);

plot3(xx_right_ql,y_q,z_bound*ones(size(y_q)),'--','LineWidth',2,'Color',c_50pro)

%---------------------
ind_X=repmat((1:81),36,1);
ind_pro=ind_X<l_bound_titer_all;
ind_brk=ind_X>r_bound_titer_all;
ind_pit=ind_X>=l_bound_titer_all|ind_X<=r_bound_titer_all;

min_pro=min(lg2_titer_2(ind_pro));
max_pro=max(lg2_titer_2(ind_pro));

min_pit=min(lg2_titer_2(ind_pit));
max_pit=max(lg2_titer_2(ind_pit));

min_brk=min(lg2_titer_2(ind_brk));
max_brk=max(lg2_titer_2(ind_brk));


ind_X=repmat(x_q,length(y_q),1);
ind_pro=ind_X<xx_left_q';
ind_brk=ind_X>xx_right_q';
ind_pit=ind_X>=xx_left_q'|ind_X<=xx_right_q';


CC=0.75*ones([size(XX),3]);


c_pro=[linspace(230,139,256)',linspace(225,189,256)',linspace(226,164,256)']/256;
c_pit=[linspace(255,216,256)',linspace(230,121,256)',linspace(218,127,256)']/256;
c_brk=[linspace(230,118,256)',linspace(256,157,256)',linspace(256,234,256)']/256;



for i=(1:length(y_q))
    for j=(1:length(x_q))
        val=lg2_titer_2(ceil(i/10),ceil(j/2));
        if ind_pro(i,j)
            norm_value=(val-min_pro)/(max_pro-min_pro);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_pro(c_ind,:);
            
        elseif ind_brk(i,j)
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

surface(ax_inten,XXl,YY,z_bound*ones(size(XX)),CC,'EdgeColor','none','FaceLighting','none','FaceAlpha',0.5)
surface(ax_inten,XXl,YY,-1*ones(size(XX)),'EdgeColor','none','FaceColor',c_50pro,'FaceAlpha',0.25)

l_bound_titer=l_bound_titer_all(10);
r_bound_titer=r_bound_titer_all(10);
y_tex_inten=11*30;
z_tex=z_bound+1;
fsize=12;
text(ax_inten,-2^(-distance(l_bound_titer)/2),y_tex_inten,z_tex,{'Immune-Imprinting- ';'Protection Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.45,0],'FontSize',fsize,'FontWeight','bold')


text(ax_inten,-2^(-(distance(l_bound_titer-3)+distance(r_bound_titer-3))/2),y_tex_inten,z_tex,{'Immune-Imprinting-';'Pitfall Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.6,0,0],'FontSize',fsize,'FontWeight','bold')


text(ax_inten,-2^(-(ad_range(2)+distance(r_bound_titer))/2),y_tex_inten,z_tex,{'Immune-Imprinting-';'Breakthrough Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.45],'FontSize',fsize,'FontWeight','bold')

%text(ax_inten,-0.0025,10*30,-0.8,'50% escape','Color',0.7*c_50pro,'VerticalAlignment','bottom','HorizontalAlignment','center')
view(ax_inten,[25.1188021427944 21.9067415730338]);