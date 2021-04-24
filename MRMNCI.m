
function [T,MRMRFS]=MRMNCI(X1,Y1)
mex computeCMImatrix_4.cpp;
mex computeMImatrix_4.cpp;
data=X1;
[nn,mm]=size(data);
rows=(1:nn);
traincount=floor((0.9)*nn);
train_rows=randsample(rows,traincount);
data1=data(train_rows,:);
C=Y1(train_rows);
data=normalize_max(data1); %normalize the data to the [-1,1] range
a=myQuantileDiscretize(data,5); %discretize the data to 5 equal-frequency bins
T=a;
T(:,mm+1)=C;
[MRMRFS]=myQPFS(a,C);
end


