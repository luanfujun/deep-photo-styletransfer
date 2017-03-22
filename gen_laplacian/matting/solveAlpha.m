function alpha=solveAlpha(I,consts_map,consts_vals,varargin)
  
  [h,w,c]=size(I);
  img_size=w*h;

 

  A=getLaplacian1(I,consts_map,varargin{:});
  
 

  D=spdiags(consts_map(:),0,img_size,img_size);
  lambda=100;
  x=(A+lambda*D)\(lambda*consts_map(:).*consts_vals(:));
 

  alpha=max(min(reshape(x,h,w),1),0);
