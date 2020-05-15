function [B] = ALLLboost(B,d)
% boosted LLL reduction in imaginary quadratic fields
% parameters: input: basis B; ring parameter: d; (e.g., d=3 for Eisenstein; d=1 for Gaussian)
% output: reduced basis B;
% author: Shanxiang Lyu, shanxianglyu@gmail.com

if nargin==1
d=1;
end
[Q,R]=qr(B);
[m,n]=size(B);
delta=0.99;
i=2;
while i<=n
    
    Ri=R(:,i);
    %size reduction
    for k=i-1:-1:1 
     q=quan(R(k,i)/R(k,k),d);    
          if q~=0  %if q is nonzero
              R(:,i)=R(:,i)-q*R(:,k);
          end    
    end

    if norm(R(:,i))>norm(Ri) && quan(Ri(i-1)/R(i-1,i-1),d)==0 %if the non-reduced Ri is shorter && effective
    R(:,i)=Ri;
    end
    
    if delta*abs(R(i-1,i-1))^2>abs(R(i,i))^2+abs(R(i-1,i))^2  %Lovasz fails
        V=[0 1;1 0];   
        %Givens matrix
        tempsum=sqrt(abs(R(i-1,i))^2+abs(R(i,i))^2);
        alpha=R(i-1,i)/tempsum;
        beta=R(i,i)/tempsum;
        R(:,i-1:i)=R(:,i-1:i)*V;
        G=eye(m);
        G(i-1:i,i-1:i)=[alpha',beta';-beta,alpha];
        R=G*R;
        Q=Q*(G)';
        i=max(i-1,2);
    else
        i=i+1;
    end   
end
B=Q*R;
end

function r= quan(x,d)
% quantizing a complex number to rings of imaginary quadratic fields
% author: Shanxiang Lyu, shanxianglyu@gmail.com

if nargin==1
    d=1;
end

if mod(-d,4)==1 %%TYPE II
    
    r1=round(real(x))+sqrt(-d)*round(imag(x)/sqrt(d));
    s=.5+.5*sqrt(-d);
 
    r2=round(real(x-.5))+sqrt(-d)*round(imag(x-.5*sqrt(-d))/sqrt(d))+s;
    
    if abs(x-r1)<abs(x-r2)
        r=r1;
    else
        r=r2;
    end
else %% TYPE I
    r=round(real(x))+sqrt(-d)*round(imag(x)/sqrt(d));
end
end