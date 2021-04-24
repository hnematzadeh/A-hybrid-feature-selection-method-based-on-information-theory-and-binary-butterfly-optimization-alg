

function [SPECCMI_Fs,MRMTRFS,QPFS_JMImat,weights,CMIMFS,CIFEFS,timeQPJMI,timeCMIM,timeCIFE,timeJMI]=mySPEC_CMI(a,C,maxFeature,QPFS_JMImat)
[n dim]=size(a);
if nargin<3 maxFeature=dim;end;

H=zeros(dim,dim);
f=zeros(dim,1);

tic
fprintf('QPFS-CMI started, calculating the MI matrix...\n')
if nargin<4||isempty(QPFS_JMImat) tic;
    H=computeCMImatrix_4([a C]); %faster version
    matrixTime=toc;QPFS_JMImat=H; 
else H=QPFS_JMImat;end;


tic;MRMTRFS=mRMTR(H,maxFeature);timeJMI=toc()+matrixTime;
tic;CMIMFS=CMIM(H,maxFeature);timeCMIM=toc()+matrixTime;
toc;CIFEFS=CIFE(H,maxFeature);timeCIFE=toc()+matrixTime;


tic;
fprintf('Solving the QPFS-JMI formulation...\n');
H=(H+H')/2;
[V,D] = eigs(H,1);

x=V(:,1);
x=abs(x/norm(x));
y=zeros(dim,2);
y(:,1)=-x;
for i=1:dim y(i,2)=i;end;
y=sortrows(y);
y=[y(:,2) -y(:,1)];
SPECCMI_Fs=y(:,1);
weights=y(:,2);
timeQPJMI=toc()+matrixTime;


end

%doing MR-MTR feature selection on a pre-computed CMI matrix H
function [best_fs]=mRMTR(H,maxFeature)
[dim dim]=size(H);
if nargin<2 maxFeature=dim;end;

fprintf('MR-MTR started...\n');
%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=H(i,i);
    if CMI>max_MI
        max_MI=CMI;
        firstFeature=i;
    end
end

for i=1:dim %create the JMI matrix
    for j=1:dim
        if j==i continue;end;
        H(i,j)=H(i,i)+H(i,j);  %I(Xi;C) + I(Xj;C|Xi)= I(Xi,Xj;C)
    end
end

best_fs=zeros(1,maxFeature);
best_fs(1)=firstFeature;
%fprintf('Adding first feature %d\n',firstFeature);

selected=zeros(1,dim);
selected(best_fs(1))=1;

for j=2:maxFeature
    max_inc=-inf;
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       totalJMI=sum(H(i,best_fs(1:j-1)));
       if totalJMI>max_inc
           max_inc=totalJMI;
           bestFeature=i;
       end
    end
    
    best_fs(j)=bestFeature;
    selected(bestFeature)=1;
%     fprintf('Adding %d\n',bestFeature);    
end
end

%performing conditional MI maximization (CMIM)
function best_fs=CMIM(H,maxFeature)
H=H'; %so that H(i,j)=I(X_i;C|Xj)??

[dim dim]=size(H);
if nargin<2 maxFeature=dim;end;

fprintf('CMIM started...\n');
%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=H(i,i);
    if CMI>max_MI
        max_MI=CMI;
        firstFeature=i;
    end
end

best_fs=zeros(1,maxFeature);
best_fs(1)=firstFeature;
%fprintf('Adding first feature %d\n',firstFeature);

selected=zeros(1,dim);
selected(best_fs(1))=1;

for j=2:maxFeature
    max_red=-inf; %max of min conditional relevancy
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       [mini]=min(H(i,best_fs(1:j-1)));
       if max_red<mini
           max_red=mini;
           bestFeature=i;
       end
    end
    
    best_fs(j)=bestFeature;
    selected(bestFeature)=1;
%     fprintf('Adding %d\n',bestFeature);    
end

end



%conditional infomax (CIFE)
function best_fs=CIFE(H,maxFeature)
H=H'; %so that H(i,j)=I(X_i;C|Xj)??

[dim dim]=size(H);
if nargin<2 maxFeature=dim;end;

fprintf('CIFE started...\n');
%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=H(i,i);
    if CMI>max_MI
        max_MI=CMI;
        firstFeature=i;
    end
end

best_fs=zeros(1,maxFeature);
best_fs(1)=firstFeature;
%fprintf('Adding first feature %d\n',firstFeature);

selected=zeros(1,dim);
selected(best_fs(1))=1;

for j=2:maxFeature
    bestobj=-inf; %minimum conditional redundancy
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       obj=(2-j)*H(i,i)+ sum(H(i,best_fs(1:j-1))); %CIFE criterion, as in PR paper
       if bestobj<obj
            bestobj=obj;
           bestFeature=i;
       end
    end
    
    best_fs(j)=bestFeature;
    selected(bestFeature)=1;
%     fprintf('Adding %d\n',bestFeature);    
end

end