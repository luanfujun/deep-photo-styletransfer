function alpha=solveAlphaC2F(I,consts_map,consts_vals,levels_num,active_levels_num,varargin);

 
  
  
  if (length(varargin)>0)
    if (~isempty(varargin{1}))
      thr_alpha=varargin{1};
    end
  end  
 
  if (~exist('thr_alpha','var'))
    thr_alpha=0.02;
  end
  erode_mask_w=1;
  active_levels_num=max(active_levels_num,1);
  if (levels_num>1)
    sI=downSmpIm(I,2);
    s_consts_map=round(downSmpIm(double(consts_map),2));
    s_consts_vals=round(downSmpIm(double(consts_vals),2));
    s_alpha=solveAlphaC2F(sI,s_consts_map,s_consts_vals,levels_num-1,...
                          min(levels_num-1,active_levels_num),varargin{:});
    alpha=upSampleAlphaUsingImg(s_alpha,sI,I,varargin{2:end});
    talpha=alpha.*(1-consts_map)+consts_vals;
    consts_map=min(consts_map+imerode((alpha>=(1-thr_alpha)),ones(erode_mask_w*2+1))+imerode((alpha<=(thr_alpha)),ones(erode_mask_w*2+1)),1);
    consts_vals=round(talpha).*consts_map;
    %figure, imshow([consts_map,alpha])
  end
  if (active_levels_num>=levels_num)
    alpha=solveAlpha(I,consts_map,consts_vals,varargin{2:end});
  end
  

