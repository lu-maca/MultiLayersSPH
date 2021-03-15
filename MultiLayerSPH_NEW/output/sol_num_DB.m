clc
clear
close all

f = @(cm) -8*9.81*0.5*cm.^2.*(sqrt(9.81*1)-cm).^2 + (cm.^2-9.81*0.5).^2.*(cm.^2+9.81*0.5);
ans = 2.6704%0.6682;
cm = ans;


tempi = [50];
t =0;
j=0;
conta = 0;
figure
for i = tempi
  t= 0.49887520459924045;
  
  filename = strcat("file",int2str(i),".dat")

  tab = dlmread(filename);
  x = tab(:,1);
  h = tab(:,2);
  
  
  plot(x,h)
  ylim([0.4 1.1])
  xlim([3 7])
  clear h1
  x1 = [0:0.001:10];
  h0 = 1;
  c = sqrt(9.8*h0);
  n=0;
  for idx = 1: length(x1)
      m = x1(idx);
      if (m>= 0 && m< 5- t*c)
          h1(idx) = h0;
      elseif (m>= 5- t*c && m<5+(2*c-3*cm)*t)
          h1(idx) = 4/(9*9.81)*(c-(m-5)/(2*t))^2;
      elseif (m>= 5+(2*c-3*cm)*t && m< 5+t*(2*cm^2*(c-cm))/(cm^2-9.81*0.5))
          h1(idx) = cm^2/9.81;
      elseif (m>= 5+t*(2*cm^2*(c-cm))/(cm^2-9.81*0.5))
          h1(idx) = 0.5;
      end
  end

  disp(t)

  hold on

  plot(x1,h1)
  plot(x1,h1)
  plot(x1,h1)
  

  %figure

  pause(0.0001)

 end
