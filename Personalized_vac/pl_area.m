function pl_area(x,y,colors,edge_alpha,axisYLim,YLab,XTick_width,Lable_off)

    ax=gca;
    ax_width=1;
    % 绘制area图
    a=area(ax,x,y','EdgeAlpha',edge_alpha);
    hold(ax,'on');
    plot(ax,[min(x),max(x)],[axisYLim(2),axisYLim(2)],'k','LineWidth',ax_width)
    plot(ax,[max(x),max(x)],[axisYLim(1),axisYLim(2)],'k','LineWidth',ax_width)
    ax.LineWidth=ax_width;
    for k=(1:length(a))
        a(k).FaceColor=colors(k,:);
    end
    set(ax, 'TickLength', [0.003, 0.003]);

    ax.XAxisLocation='bottom';
    ax.YAxisLocation='left';
    ylabel(ax,YLab)
    
    ax.XTick=(min(x):XTick_width:max(x));
    ax.YLim=axisYLim;
    ax.XLim=[min(x),max(x)];
    ax.FontSize=10;
    ax.YLabel.FontSize=12;
    set(ax,'TickDir','out')
    box(ax,'off')
    if Lable_off
        ax.XTickLabel=[];

    end
    
end
