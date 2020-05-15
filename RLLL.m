function [B2r] = RLLL(Bc,d)
% LLL reducing a complex basis (n-dim) via a transform to the real field (2n-dim)
% parameters: input: basis Br; ring parameter: d; (e.g., d=3 for Eisenstein; d=1 for Gaussian)
% output: reduced basis Bc;
% author: Shanxiang Lyu, shanxianglyu@gmail.com

if nargin==1
d=3;
end
n=size(Bc,2);
if mod(-d,4)==1
    Phi=[1 .5;0 sqrt(d)/2];
else
    Phi=[1 0;0 sqrt(d)];
end
Br=[real(Bc),-imag(Bc);imag(Bc),real(Bc)];%general complex basis
Bd=Br*(kron(Phi,eye(n)));%specific basis w.r.t. d
B2r =LLL(Bd);
