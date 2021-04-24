
function [MRMRFS]=myQPFS(a,C,maxFeature,QPFS_MImat)
[n dim]=size(a);
if nargin<3 maxFeature=dim;end;
H=zeros(dim,dim);
HH=zeros(dim,dim);
f=zeros(dim,1);
maxE=log(max(max(a)));
if nargin<4||isempty(QPFS_JMImat) tic;
    HH=computeCMImatrix_4([a C]); %faster version
    matrixTime=toc;QPFS_JMImat=HH; 
else HH=QPFS_JMImat;end;

fprintf('QPFS started, calculating the MI matrix...')

if nargin<4||isempty(QPFS_MImat) 
    tic;
    H=computeMImatrix_4([a C]);    
    matrixTime=toc;
    QPFS_MImat=H; 
else H=QPFS_MImat;end;

f=H(1:end-1,end);
H=H(1:end-1,1:end-1);

%-doing maxRel and MRMR
tic;[MRMRFS]=mrmr(H,HH,f,maxFeature);timeMRMR=toc()+ matrixTime*(maxFeature/dim ); %MRMR



end


function [MRMRFS]=mrmr(H,HH,f,maxFeature)
[N dim]=size(H);
if nargin<3 maxFeature=dim;end;

fprintf('MRMR started...\n');

%select first feature as the one with max MI with C
max_MI=0;firstFeature=1;
for i=1:dim
    CMI=f(i);
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

%forward selection
%fprintf('Forward selection...\n');
for j=2:maxFeature
    max_inc=-inf;
    bestFeature=0;
    for i=1:dim
       if selected(i) continue;end;
       rel=f(i);       
       red=sum(H(i,best_fs(1:j-1)))/(j-1);
       mini=min(HH(i,best_fs(1:j-1)));
        inc=rel-red+mini;
       if inc>max_inc
           max_inc=inc;
           bestFeature=i;
       end
    end
    
    best_fs(j)= bestFeature;
    selected(bestFeature)=1;
    %fprintf('Adding %d\n',bestFeature);    
end
MRMRFS=best_fs;
end


