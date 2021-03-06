clc
clear
close all
step=20; 
dt = 0.01;

tempi =[0:step:35000];%
t = 0-dt*step;
j=0;
conta = 0
for i = tempi
  t= t+dt*step

  filename = strcat("file",int2str(i),".dat")

    conta=conta+1;
  tab = dlmread(filename);
  x = tab(:,1);
  h = tab(:,2);
  h_pos(conta) = h(20);
  tempo(conta) = t;


   plot(x,h)
   ylim([0.00 0.2])

  pause(0.0001)

endfor
figure
plot(tempo,h_pos)

