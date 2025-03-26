function fig3

close all;

L = 5.0;            % Domain length (one dimension)
N = 20;             % Number of lines to plot

% Dirichlet problem
n_grid = 1:1:N;

dirichlet_eigenvalues = (pi*n_grid/L).^2;
c_grid = linspace(-10,10,500);

figure(1)
h = plot(c_grid, c_grid./(dirichlet_eigenvalues(:)),'linewidth',3.0);
h(1).Color = [1 0 0];
for i=2:N
    h(i).Color = [0 0 1];
end
hold on
z = xline(0,'-.','$\lambda_0$ (Neumann case only)','interpreter','latex');
z.LineWidth = 3.0;
z.LabelHorizontalAlignment = 'left';
z.LabelVerticalAlignment = 'middle';
z.FontSize = 24.0;
text(2.75,9.25,'$\lambda_1$','interpreter','latex','fontsize',24.0);
text(8.25,6,'$\lambda_2$','interpreter','latex','fontsize',24.0);
text(8.25,2.75,'$\lambda_3$','interpreter','latex','fontsize',24.0);
text(8.75,1.75,'$\lambda_4$','interpreter','latex','fontsize',24.0);
hold off

legend('$\mathcal{A}_1 = \{(c,d) : c = \lambda_1 d \}$','$\mathcal{A}_n = \{(c,d) : c = \lambda_n d \}$','location','northwest','interpreter','latex','fontsize',26.0)

% Styling stuff
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = [0.5, 0.5, 0.5]; % Custom grid color
    ax.GridLineStyle = '-';
    ax.Layer = 'top'; % Make sure the grid is on top
    
    % Remove X and Y labels while keeping the grid lines
    ax.FontSize = 14;
    ax.GridAlpha = 0.15;
    ax.LineWidth = 2;
    % Thicken tick marks
    ax.TickLength = [0.01, 0.01];  % Adjust the tick length to make them more prominent
    ax.TickDir = 'in';  % Makes tick marks point outwards for a cleaner look

ylim([ 0 c_grid(end)]);
xlim( [c_grid(1) c_grid(end)  ] );
axis square;
%

% Labels etc.
xlabel('$c = c_1 - c_2$','interpreter','latex','fontsize',36.0);
ylabel('$d = d_1 - d_2$','interpreter','latex','fontsize',36.0);

b_grid = linspace(-5,5,501);
d_grid = linspace(0.5,10,501);
c1 = zeros(length(b_grid),length(d_grid));
c2 = c1;
c3 = c1;

for i=1:length(b_grid)
    for j=1:length(d_grid)
    c1(j,i) = b_grid(i)*b_grid(i)/d_grid(j) + d_grid(j)*(1*pi/L)^2;
    c2(j,i) = b_grid(i)*b_grid(i)/d_grid(j) + d_grid(j)*(2*pi/L)^2;
    c3(j,i) = b_grid(i)*b_grid(i)/d_grid(j) + d_grid(j)*(3*pi/L)^2;
    end
end

figure(2)
h1 = surf(b_grid,d_grid,c1, 'FaceColor','r', 'FaceAlpha',0.9, 'EdgeColor','none');
hold on
h2 = surf(b_grid,d_grid,c2, 'FaceColor','y', 'FaceAlpha',0.75, 'EdgeColor','none');
h3 = surf(b_grid,d_grid,c3, 'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','none');
hold off
grid on

zl = zlabel('$c = c_1 - c_2 = \lambda_n(d,b)$','interpreter','latex','fontsize',30.0);
xl = xlabel('$b = b_1 - b_2$','interpreter','latex','fontsize',30.0);
yl = ylabel('$d = d_1 - d_2$','interpreter','latex','fontsize',30.0);

legend([h1, h2, h3], {'$c = \lambda_1(d,b)$', '$c = \lambda_2(d,b)$', '$c = \lambda_3(d,b)$'}, 'Interpreter', 'latex', 'FontSize', 24);

% Styling stuff, same as above 
ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = [0.5, 0.5, 0.5]; % Custom grid color
    ax.GridLineStyle = '-';
    ax.Layer = 'top'; % Make sure the grid is on top

ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.ZAxis.FontSize = 14; % For 3D plots
zl.FontSize = 36;
yl.FontSize = 36;
xl.FontSize = 36;
    
    ax.GridAlpha = 0.15;
    ax.LineWidth = 2;
    
    % Thicken tick marks
    ax.TickLength = [0.01, 0.01]; 
    ax.TickDir = 'in';  
    
hold off
end