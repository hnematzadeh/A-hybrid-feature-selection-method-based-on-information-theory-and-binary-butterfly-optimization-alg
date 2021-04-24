clc;
clear all;
load('CNS.mat');
dsize=size(CNS,2)-1;
X=CNS(:,1:dsize-1);
Y=CNS(:,dsize);
dsize1=floor(0.2*(dsize));
rank=zeros(5,dsize1);
rank1=zeros(1,dsize1);
rank2=zeros(1,dsize1);
rank3=zeros(1,dsize1);

    
    [T1,MRMRFS1]=MRMNCI(X,Y);
    [dim1,dim2]=size(T1);
    rank1(nk,:)=MRMRFS1(1,1:dsize1)
    for i=1 :dsize1
      ii=MRMRFS1(i);
      for j=1:dim1
         Temp7(j,i)=T1(j,ii);
      end
    end
    
    Temp7(:,dsize1+1)=T1(:,dim2);
    %save(Temp7,'Temp7');
    
    
    
    SearchAgents_no=10; % Number of search agents
    Max_iteration=100; % Maximum number of iterations
    
    [XX]=initialization(SearchAgents_no,dsize1);
    for i=1:SearchAgents_no
         fobj(i)=Get_Functions_details(XX(i,:),Temp7,dsize1);
    end
    [Best_score,Best_pos,cg_curve]=BOA(SearchAgents_no,Max_iteration,dim,fobj,XX,Temp7);
    
    
    sel=find(Best_pos==1);
    Nsel=size(sel,2)
    for k=1:Nsel
        selItem=sel(k);
        rank2(k)=rank1(selItem);
    end
    data=Temp7(:,sel);
    YY=Temp7(:,dim2);
    [fishRank,fishScore] = ftSel_fisher(data,YY);
    for i=1:Nsel
         Fisher(fishRank(i))=i;
    end
    [ranks,weights] = relieff(data,YY,10);
    for i=1:Nsel
         REL(ranks(i))=i;
    end
    Re(1,:)=Fisher(1,:);
    Re(2,:)=REL(1,:);
    for i=1:Nsel
         res(i)=(Re(1,i)/Nsel)*(Re(2,i)/Nsel);
    end
    [ResScore,ResRank] =sort(res,'ascend') ;
    for k=1:30
        TR=ResRank(k);
        rank3(k)=rank2(TR);
    end
    Rank(nk,:)=rank3(1,:);

  