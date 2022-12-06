function  [] = PlotFig(U,ExTriCnt,kind,tt)
% This function plots numerical results

% Inputs:
% U: cell average values at mesh triangles at time tt
% ExTriCnt: barycenter of triangles
% kind of method ('hybrid' or 'weno')
% tt: time level
%
if strcmp(kind,'hybrid')
    P = ExTriCnt;
    Val=U;
    P_x = linspace(min(P(:,1)), max(P(:,1)), 100);
    P_y = linspace(min(P(:,2)), max(P(:,2)), 100);
    [X,Y] = meshgrid(P_x, P_y);
    Z = griddata(P(:,1),P(:,2),Val,X,Y);
    figure('NumberTitle','on')  
    surf(X, Y, Z);
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    grid on
    set(gca,'TickLabelInterpreter','latex')
    set(gcf, 'Position', [300 300 400 400])
    set(gca, 'XTick', [-0.5 0 0.5])
    set(gca, 'YTick', [-0.5 0 0.5])
    title(num2str(sprintf('Hybrid, t = %1.1f',tt)))
    view(45,45)
 end
if strcmp(kind,'weno')
    P = ExTriCnt;
    P_x = linspace(min(P(:,1)), max(P(:,1)), 100);
    P_y = linspace(min(P(:,2)), max(P(:,2)), 100);
    [X,Y] = meshgrid(P_x, P_y);
    Z = griddata(P(:,1),P(:,2),U,X,Y);
    figure('NumberTitle','on') 
    surf(X, Y, Z);
    xlim([-0.5,0.5])
    ylim([-0.5,0.5])
    grid on
    set(gcf, 'Position', [300 300 400 400])
    set(gca, 'XTick', [-0.5 0 0.5])
    set(gca, 'YTick', [-0.5 0 0.5])
    title(num2str(sprintf('WENO, t = %1.1f',tt)))
    view(45,45)
end