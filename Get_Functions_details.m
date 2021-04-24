
function o= Get_Functions_details(XXF,Temp,d)



alpha=0.99;
beta=0.001;
gama=1-(alpha+beta);


sel=find(XXF==1);
Nsel=size(sel,2);
TT=Temp(:,sel);
YY=Temp(:,d+1);
o=alpha*accu(TT,YY)+beta*((d-Nsel)/d)+gama*computegains(TT,YY);

end



