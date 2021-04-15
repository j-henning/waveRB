% Calculates the product (A (x) B) x 
% The calculatation of the above product is equivalent to
% (B * X *  A') where a is a matrix containing the vector x
% with the desired dimensions



function b = kronVectorProduct(A, B, x)
if nargin == 3
    x_matrix = reshape(x, [ size(B,2) size(A,2)]);
    b = reshape(B* x_matrix * A', [size(B,1)*size(A,1) 1]);
end
end