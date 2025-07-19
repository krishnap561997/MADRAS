clc, clear all, close all
fileID = fopen('input.txt','r');
C = textscan(fileID,'%f %f %f %f','Delimiter','\n');

nf       = C{1}
factor   = C{2} % This is the same as the restart time in par file
ifpar    = C{3} % This is the same as last .h5 file in restart folder 
ifparout = C{4} 

if(ifpar==1)
   for n=nf:-1:2 
      movefile(strcat(sprintf('par%05d.h5\n', n)),strcat(sprintf('par%05d.h5\n', n+factor)))
   end
end
if(ifparout==1)
   for n=nf:-1:2 
      movefile(strcat(sprintf('parout%05d.h5\n', n)), strcat(sprintf('parout%05d.h5\n', n+factor)))
   end
end

