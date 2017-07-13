function K = refine_posterization(I, J)
    addpath guided_filter/
    r = 8;
    eps = 0.1^2;
    I_f = I;
    for c = 1 : 3
        I_f(:,:,c) = guidedfilter_color(I, I_f(:,:,c), r, eps);       
    end 
    K = J + I - I_f;
end 