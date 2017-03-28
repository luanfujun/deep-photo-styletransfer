addpath matting/
addpath gaimc/
N = 60;

for i = 1:N
    prefix = '../examples/input/';
    in_name = [prefix 'in' int2str(i) '.png']; 
    disp(['Working on image index = ' int2str(i)]);
    
    input = im2double(imread(in_name));
    input = reshape_img(input, 700);
    size(input)
    
    close all
    figure; imshow(input);
    
    [h w c] = size(input);
    
    disp('Compute Laplacian');
    A = getLaplacian1(input, zeros(h, w), 1e-7, 1);
 
    
    disp('Save to disk');
    n = nnz(A);
    [Ai, Aj, Aval] = find(A);
    CSC = [Ai, Aj, Aval];
    %save(['Input_Laplacian_3x3_1e-7_CSC' int2str(i) '.mat'], 'CSC');
    
    [rp ci ai] = sparse_to_csr(A);
    Ai = sort(Ai);
    Aj = ci;
    Aval = ai;
    CSR = [Ai, Aj, Aval];
    save(['Input_Laplacian_3x3_1e-7_CSR' int2str(i) '.mat'], 'CSR');
 
end 
