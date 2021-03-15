clc, clear, close

u = dlmread("upart1000.dat");
u = u(find(abs(u)>0));

x = dlmread("xpart1000.dat");
x = x(find(abs(x)>0));

plot(x,u,'o')
