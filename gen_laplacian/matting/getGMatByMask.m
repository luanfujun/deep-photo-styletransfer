function [Gx,Gy,G3,G4]=getGMatByMask(w,h,mask)
imgSize=w*h;

dS=[1,-1];
filtSizeS=1;


%indsGx1=[]; indsGx2=[]; valsGx=[];
%indsGy1=[]; indsGy2=[]; valsGy=[];
indsGx1=zeros(imgSize*2,1);
indsGx2=zeros(imgSize*2,1);
valsGx=zeros(imgSize*2,1);
indsGy1=zeros(imgSize*2,1);
indsGy2=zeros(imgSize*2,1);
valsGy=zeros(imgSize*2,1);
indsG31=zeros(imgSize*2,1);
indsG32=zeros(imgSize*2,1);
valsG3=zeros(imgSize*2,1);
indsG41=zeros(imgSize*2,1);
indsG42=zeros(imgSize*2,1);
valsG4=zeros(imgSize*2,1);

indy=0; indx=0; ind3=0; ind4=0;

for x=1:w-1,	
  for y=1:h,
    if ((~mask(y,x))&(~mask(y,x+1)))
      continue
    end  
    for disp=0:filtSizeS,
      indx=indx+1;
      indsGx1(indx)=imIndexToVect(y,x,h);
      indsGx2(indx)=imIndexToVect(y,x+disp,h);
      valsGx(indx)=dS(disp+1);
    end
  end
end
for x=1:w,	
  for y=1:h-1,
    if ((~mask(y,x))&(~mask(y+1,x)))
      continue
    end
    for disp=0:filtSizeS,
      indy=indy+1;
      indsGy1(indy)=imIndexToVect(y,x,h);
      indsGy2(indy)=imIndexToVect(y+disp,x,h);
      valsGy(indy)=dS(disp+1);
    end;
  end;
end

for x=1:w-1,	
  for y=1:h-1,
    if ((~mask(y,x))&(~mask(y+1,x+1)))
      continue
    end
    for disp=0:filtSizeS,
      ind3=ind3+1;
      indsG31(ind3)=imIndexToVect(y,x,h);
      indsG32(ind3)=imIndexToVect(y+disp,x+disp,h);
      valsG3(ind3)=dS(disp+1);
    end
  end
end

for x=1:w-1,	
  for y=2:h,
    if ((~mask(y,x))&(~mask(y-1,x+1)))
      continue
    end
    for disp=0:filtSizeS,    
      ind4=ind4+1;
      indsG41(ind4)=imIndexToVect(y,x,h);
      indsG42(ind4)=imIndexToVect(y-disp,x+disp,h);
      valsG4(ind4)=dS(disp+1);
    end;
  end;
end;
%'done inds'
indsGx1=indsGx1(1:indx);
indsGx2=indsGx2(1:indx);
valsGx=valsGx(1:indx);
indsGy1=indsGy1(1:indy);
indsGy2=indsGy2(1:indy);
valsGy=valsGy(1:indy);
indsG31=indsG31(1:ind3);
indsG32=indsG32(1:ind3);
valsG3=valsG3(1:ind3);
indsG41=indsG41(1:ind4);
indsG42=indsG42(1:ind4);
valsG4=valsG4(1:ind4);



Gx=sparse(indsGx1,indsGx2,valsGx,imgSize,imgSize);
Gy=sparse(indsGy1,indsGy2,valsGy,imgSize,imgSize);
G3=sparse(indsG31,indsG32,valsG3,imgSize,imgSize);
G4=sparse(indsG41,indsG42,valsG4,imgSize,imgSize);
