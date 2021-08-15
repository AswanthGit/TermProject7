clear all;clc
load stkadata

%% defining the lengths of the links
l1 = 8.7; l2 = 8.7;
l3 = 4; l4 = 4; l5 = 4; l6 = 4;
l7 = 3;

figure
axis([-9 9 -9 9]);
N = size(t,2);
phi1 = pcoordsall(3,1);
phi2 = pcoordsall(6,1);
phi3 = pcoordsall(9,1);
phi4 = pcoordsall(12,1);
phi5 = pcoordsall(15,1);
phi6 = pcoordsall(18,1);
phi7 = pcoordsall(21,1);

xA = -3; yA = 0;
xB = xA+l1*cos(phi1);yB = yA+l1*sin(phi1);
xC = l7*cos(phi7); yC = l7*sin(phi7);
xD = xA+l2*cos(phi2);yD = yA+l2*sin(phi2);
xE = xB+l4*cos(phi4); yE = yB+l4*sin(phi4);

line([7 7],[-6 6],'LineWidth',2,'LineStyle','--','Color','k');
link0 = line([0 xA],[0 yA],'LineWidth',1,'LineStyle','--','color','k');
link1 = line([xA xB],[yA yB],'LineWidth',3,'color','g');
link2 = line([xA xD],[yA yD],'LineWidth',3,'color','g');
link3 = line([xB xC],[yB yC],'LineWidth',3,'color','r');
link4 = line([xB xE],[yB yE],'LineWidth',3,'color','r');
link5 = line([xC xD],[yC yD],'LineWidth',3,'color','r');
link6 = line([xE xD],[yE yD],'LineWidth',3,'color','r');
link7 = line([0 xC],[0 yC],'LineWidth',3,'color','b');
Apos = rectangle('Position',[xA-0.2 yA-0.2 0.5 0.5],'Curvature',[1 1],'FaceColor','r');
Opos = rectangle('Position',[-0.2 -0.2 0.5 0.5],'Curvature',[1 1],'FaceColor','r');
timedisplay = text(1,1,num2str(t(1)));

for i = 1:10:N
    set(link1,'xdata',[xA xB],'ydata',[yA yB]);
    set(link2,'xdata',[xA xD],'ydata',[yA yD]);
    set(link3,'xdata',[xB xC],'ydata',[yB yC]);
    set(link4,'xdata',[xB xE],'ydata',[yB yE]);
    set(link5,'xdata',[xC xD],'ydata',[yC yD]);
    set(link6,'xdata',[xE xD],'ydata',[yE yD]);
    set(link7,'xdata',[0 xC],'ydata',[0 yC]);
    set(timedisplay,'Position',[-6,1],'string',num2str(t(i)));
    drawnow();
%% reinitializing
phi1 = pcoordsall(3,i);
phi2 = pcoordsall(6,i);
phi3 = pcoordsall(9,i);
phi4 = pcoordsall(12,i);
phi5 = pcoordsall(15,i);
phi6 = pcoordsall(18,i);
phi7 = pcoordsall(21,i);

xA = -3; yA = 0;
xB = xA+l1*cos(phi1);yB = yA+l1*sin(phi1);
xC = l7*cos(phi7); yC = l7*sin(phi7);
xD = xA+l2*cos(phi2);yD = yA+l2*sin(phi2);
xE = xB+l4*cos(phi4); yE = yB+l4*sin(phi4);
pause(0.001);

end