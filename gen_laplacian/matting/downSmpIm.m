function I=downSmpIm(I,filtS)

if(~exist('filtS','var'))
  filtS=1;
end
if (filtS==2)
  filt=[1,4,6,4,1]/16;
end
if (filtS==1)
  filt=[1,2,1]/4;
end

for i=1:size(I,3)
  for j=1:size(I,4)
    I(:,:,i,j)=conv2(filt,filt',I(:,:,i,j),'same');
  end
end

I=I(filtS+1:2:end-filtS,filtS+1:2:end-filtS,:,:);