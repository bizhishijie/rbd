clear
wx=1.1;wy=5;N=100000;
ii=1;

fileName=sprintf('disappear/%.1f_%d_%d_%d.mat',wx,wy,N,ii);
disp(fileName);
load(fileName);
Rc=sortRc(Rc);
for ii=1:length(Rc)
    
end