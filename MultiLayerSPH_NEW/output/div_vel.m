clc,close all,clear


t=1000
time = 0.0025*t;
beta=0.0572957
nu = 0.01
htot = 0.5
pos = 25

uu = dlmread(strcat("vel",int2str(t),".dat"));
u = zeros(size(uu,1)+1,size(uu,2));
u(2:end,:) = uu;

figure
z=dlmread(strcat("height",int2str(t),".dat"));
plot(u(:,pos),[0;z(:,pos)],'o')
zz =dlmread(strcat("zcoord",int2str(t),".dat"));
zz = [0;zz(:,pos)];

for i =1:length(zz)-1
uuu(i) = 9.81/(2*nu)*sin(beta*pi/180)*1/3*(-zz(i)^2 +(3*htot-zz(i+1))*zz(i) +(3*htot-zz(i+1))*zz(i+1));
end

hold on
plot(uuu,[z(:,pos)],'kx')
plot(uuu,[z(:,pos)],'kx')


zzzz = [0:0.001:zz(end)];
uuuu = 9.81/(2*nu)*sin(beta*pi/180).*zzzz.*(2*htot-zzzz);
plot(uuuu,zzzz)
plot(uuuu,zzzz)
plot(uuuu,zzzz)
str = sprintf("time = %s s",num2str(time));
title(str);
%legend('16 layers','u_{ave} exact','Nusselt solution','Location',"northwest")

arg=0;
u = u(:,pos); u = u(2:end);
for i = 1:length(uuu)
    arg = arg + ((uuu(i)-u(i))/uuu(i))^2;
end
errL2 = sqrt(1/length(u)*arg)

xlabel('u (m/s)')
ylabel('z (m)')
ylim([0 htot+0.5*htot])
