function draw_R_square(ax,batch_titer_lgFC,batch_titer_est_lgFC)
alpha=0.5;

scatter(ax,batch_titer_lgFC(:),batch_titer_est_lgFC(:),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor','b')
hold on

min_FC=min([batch_titer_lgFC(:),batch_titer_est_lgFC(:)],[],2,'includenan');
max_FC=max([batch_titer_lgFC(:),batch_titer_est_lgFC(:)],[],2,'includenan');

min_FC=min(min_FC);
max_FC=max(max_FC);


lgFC_width=[min_FC-0.5,max_FC+0.5];

plot(ax,lgFC_width,lgFC_width,'--r','LineWidth',2)
set(ax,'XTick',(floor(min_FC):ceil(max_FC)),'YTick',(floor(min_FC):ceil(max_FC)))

mean_lgFC=mean(batch_titer_lgFC(:),'omitnan');
TSS=sum((batch_titer_lgFC(:)-mean_lgFC).^2,'omitnan');
ESS=sum((batch_titer_lgFC(:)-batch_titer_est_lgFC(:)).^2,'omitnan');
R=1-ESS/TSS;

n_data=length(batch_titer_est_lgFC(:))-sum(isnan(batch_titer_lgFC(:)));
fsize=8;
text(ax,0.05,0.82,['R^2=',num2str(round(R,4))],'Units','normalized','HorizontalAlignment','left','FontSize',fsize)
text(ax,0.05,0.9,['n=',num2str(n_data)],'Units','normalized','HorizontalAlignment','left','FontSize',fsize)
xlabel(ax,'log_2(titer FC) (Experiment) ')
ylabel(ax,'log_2(titer FC) (Model) ')
xlim(ax,lgFC_width)
ylim(ax,lgFC_width)
end