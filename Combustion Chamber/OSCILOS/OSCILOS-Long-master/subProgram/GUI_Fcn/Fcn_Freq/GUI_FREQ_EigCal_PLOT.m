
%-------------------------------------------------------------------------
function GUI_FREQ_EigCal_PLOT(varargin)
hObject             = varargin{1};
handles             = guidata(hObject);
global CI
global FDF
hAxes1              = handles.axes1;
hAxes2              = handles.axes2;
fontSize1           = handles.FontSize(1);
fontSize2           = handles.FontSize(2);
CI.EIG.pop_numMode  = get(handles.pop_numMode,  'Value');
CI.EIG.pop_PlotType = get(handles.pop_PlotType, 'Value');
if CI.EIG.APP_style~=3
    ValueSlider         = get(handles.slider_uRatio,'Value');  
else
    ValueSlider         = get(handles.slider_pRatio,'Value');
end
indexShow           = round(ValueSlider);
Eigenvalue          = CI.EIG.Scan.EigValCol{indexShow};
ValueContour        = CI.EIG.Cont.ValCol{indexShow};
pannelsize          = get(handles.uipanel_Axes,'position');
pW = pannelsize(3);
pH = pannelsize(4); 
%
guidata(hObject, handles)
%
switch CI.EIG.pop_PlotType
% {'Map of eigenvalues';
%  'Modeshape';
%  'Evolution of eigenvalue with velocity ratio'}
    case 1      % Map of eigenvalues;
    set(handles.pop_numMode,'enable','off'); 
    set(hAxes1,'position',[pW*1.5/10 pH*1.5/10 pW*7/10 pH*7/10]);  
    position_hAxes1=get(hAxes1,'position');
    try
        cbh = findobj( 0, 'tag', 'Colorbar' );
        delete(cbh)
    catch
    end
    cla(hAxes1,'reset')
    axes(hAxes1)
    hold on
    %
    contourf(hAxes1,CI.EIG.Cont.GRSp./100,CI.EIG.Cont.FreqSp,20*log10(abs(ValueContour')))
    drawnow
    ylimitUD = [CI.EIG.Scan.FreqMin CI.EIG.Scan.FreqMax];
    xlimitUD = [CI.EIG.Scan.GRMin   CI.EIG.Scan.GRMax]./100;
    hold off
    %
    set(hAxes1,'YColor','k','Box','on');
    set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    xlabel(hAxes1,  '$ Re(s)/100: \textrm{Growth rate}~~/100~~$ [rad s$^{-1}$] ',...
        'Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes1,'$ Im(s)/2\pi: \textrm{Frequency}~~$ [Hz]','Color','k',...
        'Interpreter','LaTex','FontSize',fontSize1);
    set(hAxes1,'ylim',ylimitUD,'YAxisLocation','left','Color','w');
    set(hAxes1,'xlim',xlimitUD);
    grid on
    colorbar 
    colormap(hot);
    hcb=colorbar;
    set(hcb,'Fontsize',fontSize2,'box','on','Unit','points')
    set(hcb,'position',[position_hAxes1(1)+position_hAxes1(3),...
                        position_hAxes1(2),...
                        position_hAxes1(3)./20,...
                        position_hAxes1(4).*1]);
    set(hAxes1,'position',position_hAxes1)
        hcb_ylim=get(hcb,'ylim');
        handles.hColorbar.ylimit=[  min(min(20*log10(abs(ValueContour')))),...
                                    max(max(20*log10(abs(ValueContour'))))];
        guidata(hObject, handles)
    %------------------------------------
    cla(hAxes2,'reset')
    axes(hAxes2)
    set(hAxes2,'position',get(hAxes1,'position'));
    hold on
    plot(hAxes2,real(Eigenvalue)./100,imag(Eigenvalue)./2./pi,'p',...
        'markersize',8,'color','k','markerfacecolor',[1,1,1])
    drawnow
    hold off
    set(hAxes2,     'ylim', get(hAxes1,'ylim'),...
                    'yTick',get(hAxes1,'ytick'),...
                    'yticklabel',[],...
                    'YAxisLocation','left','Color','none');
    set(hAxes2,     'xlim', get(hAxes1,'xlim'),...
                    'xTick',get(hAxes1,'xtick'),...
                    'xticklabel',[],...
                    'xcolor','b','ycolor','b','gridlinestyle','-.');
    set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',0.5)
    grid on
    %------------------------
    xt_pos=min((get(hAxes2,'xlim')))+1.15*(max((get(hAxes1,'xlim'))) - min((get(hAxes1,'xlim'))));
    yt_pos=mean(get(hAxes2,'ylim'));
    hTitle = title(hAxes2, 'Eigenvalues are located at minima');
    set(hTitle, 'interpreter','latex', 'fontunits','points','fontsize',fontSize1)
    guidata(hObject, handles);
    %--------------------
    case 2 %'Modeshape'
    set(handles.pop_numMode,'enable','on');
    set(hAxes1,'position',[pW*2.0/10 pH*1.5/10 pW*7/10 pH*3.5/10]); 
    set(hAxes2,'position',[pW*2.0/10 pH*5.0/10 pW*7/10 pH*3.5/10]); 
    try
        cbh = findobj( 0, 'tag', 'Colorbar' );
        delete(cbh)
    catch
    end
    s_star = Eigenvalue(CI.EIG.pop_numMode);             % eigenvalue
    
    switch CI.EIG.APP_style
    case {11,12}  
        [x_resample,p,u] = Fcn_calculation_eigenmode_Linear(s_star);
    case {21,22}                             % nonlinear flame model
        global HP
        HP = CI.FM.HP{CI.FM.indexMainHPinHp};
        assignin('base','HP',HP);
        FDF.num     = CI.EIG.FDF.num{indexShow};
        FDF.den     = CI.EIG.FDF.den{indexShow};
        FDF.tauf    = CI.EIG.FDF.tauf(indexShow);
        FDF.uRatio  = CI.EIG.FDF.uRatioSp(indexShow);
        assignin('base','FDF',FDF);
        [x_resample,p,u] = Fcn_calculation_eigenmode_frozen_nonlinear(s_star);
    case {3}
        CI.EIG.PR.A1minus_this_step=CI.EIG.PR.A1minusSp(indexShow);
        assignin('base','CI',CI);
        [x_resample,p,u] = Fcn_calculation_eigenmode_nonlinear_dampers(s_star);    
    end
    cla(hAxes1,'reset')
    axes(hAxes1)
    drawnow
    hold on
    for k=1:length(CI.CD.x_sample)-1
        plot(hAxes1,x_resample(k,:),abs(p(k,:)),'-','color','k','Linewidth',2)
    end
    valueLevel = round(log10(max(max(abs(p)))));
    ymax1=ceil(max(max(abs(p)))./10^valueLevel).*10^valueLevel;
    ymin1=floor(min(min(abs(p)))./10^valueLevel).*10^valueLevel;
    ylimitUD=[ymin1 ymax1+0.25*(ymax1-ymin1)];
    ytickUD=linspace(ylimitUD(1),ylimitUD(2),6);
    for ss=1:length(ytickUD)
        yticklabelUD{ss}=num2str(ytickUD(ss));
    end
    yticklabelUD{end}='';
    xmax1=max(max(x_resample));
    xmin1=min(min(x_resample));
    xlimitUD=[xmin1 xmax1];
    xtickUD=linspace(xlimitUD(1),xlimitUD(2),6);

    set(hAxes1,'YColor','k','Box','on');
    set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    xlabel(hAxes1,'$x $ [m]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes1,'$|~\hat{p}~|$ [Pa] ','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    set(hAxes1,     'ylim', ylimitUD,...
                    'yTick',ytickUD,...
                    'yticklabel',yticklabelUD,...
                    'YAxisLocation','left');
    set(hAxes1,'xlim',xlimitUD,'xtick',xtickUD);
    ylimit=get(hAxes1,'ylim');
    NSp = length(CI.CD.x_sample);
    if NSp < 10
        for k=1:length(CI.CD.x_sample)
           plot(hAxes1,[CI.CD.x_sample(k),CI.CD.x_sample(k)],ylimit,'--','linewidth',0.5,'color','k') 
        end
    end
    grid on
    hold off
    %-----------------------
    cla(hAxes2,'reset')
    axes(hAxes2)
    drawnow
    hold on
    for k=1:length(CI.CD.x_sample)-1
        plot(hAxes2,x_resample(k,:),abs(u(k,:)),'-','color','k','Linewidth',2)
    end
    set(hAxes2,'YColor','k','Box','on');
    set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    set(hAxes2,     'xlim', get(hAxes1,'xlim'),...
                    'xTick',get(hAxes1,'xtick'),...
                    'xticklabel',[]);
    xlabel(hAxes2,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes2,'$|~\hat{u}~|$ [m/s] ','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylimit=get(hAxes2,'ylim');
    if NSp < 10
    for k=1:length(CI.CD.x_sample)
       plot(hAxes2,[CI.CD.x_sample(k),CI.CD.x_sample(k)],ylimit,'--','linewidth',1,'color','k') 
    end
    end
    grid on
    hold off
    guidata(hObject, handles);
    %--------------------------------
    case 3 %'Evolution of eigenvalue with velocity ratio'
    
    if CI.EIG.APP_style==11 || CI.EIG.APP_style==12 
       set(hObject, 'String', 0);
       errordlg('There is not any evolution result','Error');
    elseif CI.EIG.APP_style==21 || CI.EIG.APP_style==22
        RatioNum = length(CI.EIG.FDF.uRatioSp);
    else
        RatioNum = length(CI.EIG.FDF.pRatioSp);    
    end
      
    Modes_Num     =length(CI.EIG.Scan.EigValCol{1});
    EigFreq       =abs(imag(CI.EIG.Scan.EigValCol{1}))/2/pi;
    frequency_gap =min(abs(diff(EigFreq)));
        if RatioNum==1
            set(hObject, 'String', 0);
            errordlg('There is not any evolution result','Error');
        else
            for M_Num=1:1:Modes_Num 
                EigFreq_mode(M_Num,1)    = abs(imag(CI.EIG.Scan.EigValCol{1}(M_Num)))/2/pi;
                EigGR_mode(M_Num,1)      = real(CI.EIG.Scan.EigValCol{1}(M_Num));
                k=2;
                Ratio_num_step=1;
                while k<=RatioNum
                      mod_fre_diff  =abs(EigFreq(M_Num)-abs(imag(CI.EIG.Scan.EigValCol{k}))/2/pi);
                      [min_fre_diff,min_mod_num]= min(mod_fre_diff);
                      if min_fre_diff<frequency_gap/2
                         EigFreq_mode(M_Num,k)     = abs(imag(CI.EIG.Scan.EigValCol{k}(min_mod_num)))/2/pi;
                         EigGR_mode(M_Num,k)       = real(CI.EIG.Scan.EigValCol{k}(min_mod_num));
                         Ratio_num_step            = Ratio_num_step+1;
                         k                         = k+1;
                      else
                          k=k+1;
                      end
                end
                Ratio_num(M_Num)                       =Ratio_num_step;
                if CI.EIG.APP_style==21 || CI.EIG.APP_style==22
                    Ratio_mode(M_Num,1:1:Ratio_num(M_Num))=CI.EIG.FDF.uRatioSp(1:1:Ratio_num(M_Num));
                else
                    Ratio_mode(M_Num,1:1:Ratio_num(M_Num))=CI.EIG.FDF.pRatioSp(1:1:Ratio_num(M_Num));
                end
                   
            end   
        end
    % Get the pop_numMode and plot its evolution 
    set(handles.pop_numMode,'enable','on');
    set(hAxes1,'position',[pW*2.0/10 pH*1.5/10 pW*7/10 pH*3.5/10]); 
    set(hAxes2,'position',[pW*2.0/10 pH*5.0/10 pW*7/10 pH*3.5/10]); 
    try
        cbh = findobj( 0, 'tag', 'Colorbar' );
        delete(cbh)
    catch
    end    
    % Calculate the eigen values for a specific mode
    x       = Ratio_mode(CI.EIG.pop_numMode,1:Ratio_num(CI.EIG.pop_numMode));
    EigFreq = EigFreq_mode(CI.EIG.pop_numMode,1:Ratio_num(CI.EIG.pop_numMode));
    EigGR   = EigGR_mode(CI.EIG.pop_numMode,1:Ratio_num(CI.EIG.pop_numMode));
    cla(hAxes1,'reset')
    axes(hAxes1)
    hold on
    plot(hAxes1,x,EigFreq,'-o','color','k','Linewidth',1)
    hold off
    ymax1=ceil(max(max(EigFreq)));
    ymin1=floor(min(min(EigFreq)));
    ylimitUD=[ymin1-0.25*(ymax1-ymin1) ymax1+0.25*(ymax1-ymin1)];
    ytickUD=linspace(ylimitUD(1),ylimitUD(2),7);
    for ss=1:length(ytickUD)
        yticklabelUD{ss}=num2str(ytickUD(ss));
    end
    yticklabelUD{end}='';
    if length(x)==1
        xmax1=max(max(x))+0.1;
        xmin1=min(min(x))-0.1;
    else
         xmax1=max(max(x));
         xmin1=min(min(x));
    end
    xlimitUD=[xmin1 xmax1];
    xtickUD=linspace(xlimitUD(1),xlimitUD(2),6);
    set(hAxes1,'YColor','k','Box','on');
    set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1);
    if CI.EIG.APP_style==11 || CI.EIG.APP_style==12 
      xlabel(hAxes1,'$\hat{u}_1/\bar{u}_1 $ [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    else
      xlabel(hAxes1,'$\hat{p}_1/\bar{p}_1 $ [$\times 10^{-5}$]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    end
    ylabel(hAxes1,'Frequency [Hz] ','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    set(hAxes1,     'ylim', ylimitUD,...
                    'yTick',ytickUD,...
                    'yticklabel',yticklabelUD,...
                    'YAxisLocation','left');
    set(hAxes1,'xlim',xlimitUD,'xtick',xtickUD);
    grid on
    hold off
    %-----------------------
    cla(hAxes2,'reset')
    axes(hAxes2)
    drawnow
    hold on
    plot(hAxes2,x,EigGR,'-o','color','k','Linewidth',1)
    hold off
    set(hAxes2,'YColor','k','Box','on');
    set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    set(hAxes2,     'xlim', get(hAxes1,'xlim'),...
                    'xTick',get(hAxes1,'xtick'),...
                    'xticklabel',[]);
    xlabel(hAxes2,'','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes2,'Growth rate [rad/s] ','Color','k','Interpreter','LaTex','FontSize',fontSize1);
    grid on
    guidata(hObject, handles);
end
assignin('base','CI',CI);                   % save the current information to the workspace


%
% -------------------------end---------------------------------------------    