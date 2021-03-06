clc,close,clear

t =300

hlayers = dlmread(strcat("zcoord",int2str(t),".dat"));
tab = dlmread(strcat("file",int2str(t),".dat"));

for i = 1:10
  plot(tab(:,1),hlayers(i,:),"Linewidth",1.5)
  hold on
end
