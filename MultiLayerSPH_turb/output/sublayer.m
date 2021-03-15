clc,close,clear

t = 200

hlayers = dlmread(strcat("zcoord",int2str(t),".dat"));
tab = dlmread(strcat("file",int2str(t),".dat"));

for i = 1:50
  plot(tab(:,1),hlayers(i,:),"Linewidth",1.5)
  hold on
end
