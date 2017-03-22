function [F,B]=solveFB(I,alpha)
  [h,w,c]=size(I);
  mask=(alpha>=0.02).*(alpha<=0.98);
  [Gx,Gy,Gd1,Gd2]=getGMatByMask(w,h,mask);
  G=[Gx;Gy;Gd1;Gd2];
  Ga=G*alpha(:);
  wgf=abs(Ga).^0.5+0.003*repmat((1-alpha(:)),4,1);
  wgb=abs(Ga).^0.5+0.003*repmat(alpha(:),4,1);

  
  wf=(alpha(:)>0.98)*100+0.03*alpha(:).*(alpha(:)>0.7)+0.01*(alpha(:)<0.02);
  wb=(alpha(:)<0.02)*100+0.03*(1-alpha(:)).*(alpha(:)<0.3)+0.01*(alpha(:)>0.98);
  
  
  wgf=spdiags(wgf(:),0,length(wgf),length(wgf));
  wgb=spdiags(wgb(:),0,length(wgb),length(wgb)); 
  wf=spdiags(wf(:),0,length(wf),length(wf));
  wb=spdiags(wb(:),0,length(wb),length(wb)); 
  
  
  for t=1:c
    tI=I(:,:,t);
    Ag=[wgf*G,sparse(size(G,1),size(G,2));sparse(size(G,1),size(G,2)),wgb* ...
        G];
    bg=zeros(size(Ag,1),1);
    Ai=[wf,sparse(w*h,w*h);sparse(w*h,w*h),wb];
    bi=[wf*tI(:).*(alpha(:)>0.02);wb*tI(:).*(alpha(:)<0.98)];
    As=[spdiags(alpha(:),0,w*h,w*h),spdiags(1-alpha(:),0,w*h,w*h)];
    bs=tI(:);
    A=[Ai;As;Ag];
    b=[bi;bs;bg];
    x=(A'*A)\(A'*b);
    F(:,:,t)=reshape(x(1:w*h),h,w);
    B(:,:,t)=reshape(x(w*h+1:end),h,w);

  end  