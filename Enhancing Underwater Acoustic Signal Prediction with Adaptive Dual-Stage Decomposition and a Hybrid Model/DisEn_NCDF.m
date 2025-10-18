function [Out_DisEn, npdf]=DisEn_NCDF(x,m,nc,MA,tau)
%
% This function calculates dispersion entropy (DisEn) of a univariate
% signal x, using normal cumulative distribution function (NCDF)
%
% Inputs:
%
% x: univariate signal - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% nc: number of classes (it is usually equal to a number between 3 and 9 - we used c=6 in our studies)
% tau: time lag (it is usually equal to 1)
%
% Outputs:
%
% Out_DisEn: scalar quantity - the DisEn of x
% npdf: a vector of length nc^m, showing the normalized number of disersion patterns of x
%
% Ref:
%
% [1] H. Azami, M. Rostaghi, D. Abasolo, and J. Escudero, "Refined Composite Multiscale Dispersion Entropy and its Application to Biomedical
% Signals", IEEE Transactions on Biomedical Engineering, 2017.
% [2] M. Rostaghi and H. Azami, "Dispersion Entropy: A Measure for Time-Series Analysis", IEEE Signal Processing Letters. vol. 23, n. 5, pp. 610-614, 2016.
%
% If you use the code, please make sure that you cite references [1] and [2].
%
% Hamed Azami, Mostafa Rostaghi, and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk, rostaghi@yahoo.com, and javier.escudero@ed.ac.uk
%
%  20-January-2017
%%


N=length(x);
sigma_x=std(x);
mu_x=mean(x);

  %% Mapping approaches

switch MA
    case   'LM'
        y=mapminmax(x,0,1);
        y(y==1)=1-1e-10;
        y(y==0)=1e-10;
        z=round(y*nc+0.5);
        
    case 'NCDF'
        y=normcdf(x,mu_x,sigma_x);
        y=mapminmax(y,0,1);
        y(y==1)=1-1e-10;
        y(y==0)=1e-10;
        z=round(y*nc+0.5);
        
    case 'LOGSIG'
        y=logsig((x-mu_x)/sigma_x);
        y=mapminmax(y,0,1);
        y(y==1)=1-1e-10;
        y(y==0)=1e-10;
        z=round(y*nc+0.5);
        
    case 'TANSIG'
        y=tansig((x-mu_x)/sigma_x)+1;
        y=mapminmax(y,0,1);
        y(y==1)=1-1e-10;
        y(y==0)=1e-10;
        z=round(y*nc+0.5);
        
    case 'SORT'
        x=x(1:nc*floor(N/nc));
        N=length(x);
        [sx osx]=sort(x);
        Fl_NC=N/nc;
        cx=[];
        for i=1:nc
            cx=[cx i*ones(1,Fl_NC)];
        end
        for i=1:N
            z(i)=cx(osx==i);
        end
        
end
        

all_patterns=[1:nc]';

for f=2:m
    temp=all_patterns;
    all_patterns=[];
    j=1;
    for w=1:nc
        [a,b]=size(temp);
        all_patterns(j:j+a-1,:)=[temp,w*ones(a,1)];
        j=j+a;
    end
end

for i=1:nc^m
    key(i)=0;
    for ii=1:m
        key(i)=key(i)*10+all_patterns(i,ii);
    end
end


embd2=zeros(N-(m-1)*tau,1);
for i = 1:m
    embd2=[z(1+(i-1)*tau:N-(m-i)*tau)]'*10^(m-i)+embd2;
end


pdf=zeros(1,nc^m);

for id=1:nc^m
    [R,C]=find(embd2==key(id));
    pdf(id)=length(R);
end

npdf=pdf/(N-(m-1)*tau);
p=npdf(npdf~=0);
Out_DisEn = -(sum(p .* log(p)))/log(nc^m);