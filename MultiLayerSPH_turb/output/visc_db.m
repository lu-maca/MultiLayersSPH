clc
clear
close all

tempi = [0:1000:100000];
t = -0.1;
j=0;
conta = 0
my_col = 0;
for i = tempi
  t= t+i*0.1*i;

  filename = strcat("file",int2str(i),".dat")

  tab = dlmread(filename);
  x = tab(:,1);
  h = tab(:,2);



  plot(x/0.44,h/0.055)

  pause(0.0001)

endfor
