function out = reshape_img(in, len)
    if nargin < 2
        len = 512;
    end 
    
    [h,w,~] = size(in);
    if h > w
        h2 = len;
        w2 = floor(w * h2 / h);
    else 
        w2 = len;
        h2 = floor(h * w2 / w);
    end 
    
    out = imresize(in, [h2 w2]);
end 