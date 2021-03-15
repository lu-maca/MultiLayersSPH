clc, clear, close all

t = [ 1:1:10 ];
c=0;
tempi = t/0.001;
for i = tempi
  filename = strcat("file",int2str(i),".dat");
  c = c+1;
  tab = dlmread(filename);
  x = tab(:,1);
  h = tab(:,2);
  for j = 1:length(h)
    if h(j) <= 1e-5
      xf(c) = x(j)*1000 -440;
      break
    end
  end
end
xf
loglog(t*((12*130*0.44^2)/(1418*9.8*0.055^3))^(-1),xf/440,'o')
grid minor
