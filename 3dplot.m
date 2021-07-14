clear 
clc
load c:/tmp/errodata3d.mat

z1=data(:,1);
i=data(:,2);
erro=data(:,3);
figure(2)
hold on
for x=1:19
    plot3(i(91*(x-1)+1:91*x),z1(91*(x-1)+1:91*x),erro(91*(x-1)+1:91*x))
end
 Legend=cell(x,1)
 for iter=1:x
   Legend{iter}=strcat('Z1=', num2str(z1(91*iter)));
 end
 legend(Legend)

hold off
xlabel('i (Z2/Z1)', 'fontweight', 'bold')
ylabel("Z1(Number of pinion's teeth", 'fontweight', 'bold')
zlabel('Difference (‰)', 'fontweight', 'bold')

axis tight
grid on
xv = linspace(min(z1), max(z1), 200);
yv = linspace(min(i), max(i), 200);
[Z1,I] = meshgrid(xv, yv);
ERRO = griddata(z1,i,erro,Z1,I);
ZERO=zeros(200,200);
figure(1) 
hold on

surf(I, Z1, ERRO);
surf(I,Z1,ZERO);

% Find the difference field.
f3 = ERRO - ZERO;
% Interpolate the difference field on the explicitly defined surface.
% Find the contour where the difference (on the surface) is zero.
contour(I, Z1, f3, [0 0],'k');


hold off
%colormap(jet(256)) 
colorbar
%Setting the axes label, type, position, rotation, etc.
xlabel('i (Z2/Z1)', 'fontweight', 'bold')
ylabel("Z1(Number of pinion's teeth", 'fontweight', 'bold')
zlabel('Difference (‰)', 'fontweight', 'bold')



yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')


zh = get(gca,'ZLabel'); % Handle of the z label
set(zh, 'Units', 'Normalized')


grid on
set(gca, 'ZLim',[-50 20])
axis tight
%camlight
lighting phong
shading interp
%set(s1,'edgecolor',[0 0 0.4],'meshstyle','both','linewidth',.15);


