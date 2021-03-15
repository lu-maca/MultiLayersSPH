clc,close all,clear

for t = [1000:10:1000]
  t=100
  time = 0.001*t;
  beta=1
  nu = 0.1
  htot = 1

  uu = dlmread(strcat("vel",int2str(t),".dat"));
  u = zeros(size(uu,1)+1,size(uu,2));
  u(2:end,:) = uu;

  tab = dlmread(strcat("file",int2str(t),".dat"));
  x = tab(:,1);
  dx = x(2)-x(1);

  for i = 1: size(u,1)
    for j = 1:size(u,2)
      if (j==1)
        dudx(i,j) = (u(i,j+1)-u(i,j))/dx;
      elseif (j==size(uu,2))
        dudx(i,j) = (u(i,j)-u(i,j-1))/dx;
      else
        dudx(i,j) = 0.5*(u(i,j+1)-u(i,j-1))/dx;
      endif
    endfor
  endfor

  hh = dlmread(strcat("height",int2str(t),".dat"));
  h = zeros(size(hh,1)+1,size(hh,2));
  h(2:end,:) = hh;

  w = zeros(size(hh,1)+1,size(hh,2));
  for i = 1: size(u,1)-1
    for j = 1:size(u,2)
      w(i+1,j) = -dudx(i+1,j)*(h(i+1,j)-h(i,j))+w(i,j);
    endfor
  endfor

  xx = zeros(size(u));
  for i =1:size(xx,1)
     xx(i,:) = x;
  endfor
  vf = quiver (xx, h, u, w)
  hold on
  plot(x,tab(:,2),"Linewidth",2)
  xlim([0 x(end)])
  ylim([0 1.2])

  figure
  z=dlmread(strcat("height",int2str(t),".dat"));
  plot(u(:,80),[0;z(:,80)],'o')
  zz =dlmread(strcat("zcoord",int2str(t),".dat"));
  zz = [0;zz(:,80)];
  %uu = 9.81/(2*0.01)*(2*0.1-zz).*zz*sin(beta*pi/180);

  for i =1:length(zz)-1
    uuu(i) = 9.81/(2*nu)*sin(beta*pi/180)*1/3*(-zz(i)^2 +(3*htot-zz(i+1))*zz(i) +(3*htot-zz(i+1))*zz(i+1));
  endfor
  red = 0.1; green = 0.; blue = 0.5;
  hold on,plot(uuu,[z(:,80)],'x','color',[red green blue]),hold on
  zzzz = [0:0.001:zz(end)];
  uuuu = 9.81/(2*nu)*sin(beta*pi/180).*zzzz.*(2*htot-zzzz);
  plot(uuuu,zzzz,"Linewidth",1.5)
  str = sprintf("time = %s",num2str(time));
  title(str);
  legend('10 strati','esatta vmedia','Location',"northwest")

  zz =dlmread(strcat("zcoord",int2str(t),".dat"));

  figure
  for i = 1:5
    hold on
    plot(x,zz(i,:),"Linewidth",1.5)
  endfor

endfor
