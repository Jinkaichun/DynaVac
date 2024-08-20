function dynaVac_pl(t,Agt,Abt,Ft,Mt,vaccine_varient,vaccine_time,varient_color,position,Yaxis_on)

    l_comb=3;
    l_pic=15;
    jj=1;
    l_sum=l_comb+4*jj+4*l_pic;
    
    %position_comb=[position(1),position(2)+position(4)*(l_sum-l_comb)/l_sum,position(3),position(4)*l_comb/l_sum];
    
    position_Ag=[position(1),position(2)+position(4)*(3*jj+3*l_pic)/l_sum,position(3),position(4)*l_pic/l_sum];
    
    position_F=[position(1),position(2)+position(4)*(2*jj+2*l_pic)/l_sum,position(3),position(4)*l_pic/l_sum];
    
    position_M=[position(1),position(2)+position(4)*(1*jj+1*l_pic)/l_sum,position(3),position(4)*l_pic/l_sum];
    
    position_Ab=[position(1),position(2)+position(4)*(0*jj+0*l_pic)/l_sum,position(3),position(4)*l_pic/l_sum];
    
    
    
    
    sub_varient=unique(vaccine_varient);
    sub_varient=sub_varient(sub_varient~=0);
    sub_varient_color=varient_color(sub_varient,:);
    edge_alpha=0.1;
    %figure(batch)
%-----------draw timeline-----------
    XTick_width=30;
    
%     ax_comb=axes;
%     set(ax_comb,'Position',position_comb)
%     hold on
%     pl_vaccine_combo_cut(t,vaccine_time,vaccine_varient,vaccine_type,varient_color,Xtrunc,XTick_width,1)
    
    
    ax_Ag=axes;
    set(ax_Ag,'Position',position_Ag)
    
    XLab='';
    if Yaxis_on
        YLab='Antigen';
    else
        YLab='';
    end
    

    axisYLim=[0.01,1.1*max(sum(Agt))];

    pl_area(t,Agt,sub_varient_color,edge_alpha,axisYLim,YLab,XTick_width,1)
    
    
    ax_Ab=axes;
    set(ax_Ab,'Position',position_Ab)
    
    XLab='';
    if Yaxis_on
        YLab='Antibody';
    else
        YLab='';
    end
    
    
    %axisYLim=Ab_lim;
    axisYLim=[0.01,1.1*max(sum(Abt))];

    pl_area(t,Abt,sub_varient_color,edge_alpha,axisYLim,YLab,XTick_width,0)

    ax_F=axes;
    set(ax_F,'Position',position_F)
    
    XLab='';
    if Yaxis_on
        YLab='Naive B cell maturity';
    else
        YLab='';
    end
    axisYLim=[0,1.1*max(sum(Ft))];
    pl_area(t,Ft,sub_varient_color,edge_alpha,axisYLim,YLab,XTick_width,1)
    
    ax_M=axes;
    set(ax_M,'Position',position_M)
    
    XLab='Time (days)';
    if Yaxis_on
        YLab='Memory B cell';
    else
        YLab='';
    end    
    axisYLim=[0,1];
    pl_area(t,Mt,sub_varient_color,edge_alpha,axisYLim,YLab,XTick_width,1)

    
end