function nI= upSmpIm(I,new_imgSize,filtS)

if(~exist('filtS','var'))
  filtS=1;
end
if (filtS==2)
  filt=[1,4,6,4,1]/8;
end
if (filtS==1)
  filt=[1,2,1]/2;
end
if (filtS==0)
  filt=[1];
end



id=floor((new_imgSize(1)-size(I,1)*2+1)/2);
iu=ceil((new_imgSize(1)-size(I,1)*2+1)/2);
jd=floor((new_imgSize(2)-size(I,2)*2+1)/2);
ju=ceil((new_imgSize(2)-size(I,2)*2+1)/2);

nI=zeros(new_imgSize(1)+2*filtS,new_imgSize(2)+2*filtS,size(I,3),size(I, ...
                                                  4));
nI(id+filtS+1:2:end-iu-filtS,jd+filtS+1:2:end-ju-filtS,:,:)=I;
nI(id+filtS-1:-2:1,:,:,:)=repmat(nI(id+filtS+1,:,:,:),ceil((id+filtS-1)/2),1);

nI(end-iu-filtS+2:2:end,:,:,:)=repmat(nI(end-iu-filtS,:,:,:),ceil((iu+filtS-1)/2),1);
nI(:,jd+filtS-1:-2:1,:,:)=repmat(nI(:,jd+filtS+1,:,:),1,ceil((jd+filtS-1)/2));
nI(:,end-ju-filtS+2:2:end,:,:)=repmat(nI(:,end-ju-filtS,:,:),1,ceil((ju+filtS-1)/2));




for i=1:size(nI,3)
  for j=1:size(nI,4)
    nI(:,:,i,j)=conv2(filt,filt',nI(:,:,i,j),'same');
  end
end



nI=nI(filtS+1:end-filtS,filtS+1:end-filtS,:,:);



