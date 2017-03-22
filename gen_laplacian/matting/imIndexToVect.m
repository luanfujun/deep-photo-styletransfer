% function I = imIndexToVect(Y,X,imHeight)
function I = imIndexToVect(Y,X,imHeight)

I = reshape(Y + (X-1)*imHeight,prod(size(X)),1); 