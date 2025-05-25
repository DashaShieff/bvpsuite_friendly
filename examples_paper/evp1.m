function [ret]=evp1(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

%Inputs for bvpsuite
%Programmed by Georg Kitzhofer, July 2004

%Kontrollnummer43753976430976655145
if (nargin<15) zc=[]; if (nargin<14) DpRpar=[]; if (nargin<13) DpRgl=[]; if (nargin<12) Dppar=[]; if (nargin<11) Dpgl=[]; if (nargin<10) p=[]; if (nargin<9) D2=[]; if (nargin<8) D1=[]; if (nargin<7) zb=[]; if (nargin<6) za=[]; if (nargin<5) z=[]; if (nargin<4) t=[]; if (nargin<3) R=[];
if (nargin<2) D=[];end;end;end;end;end;end;end;end;end;end;end;end;end;end;
ret=intern(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc);

function ret=intern(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)
switch nr
     case 'x0'
         for j=1:(length(intern('x1'))-1)*(intern('n')*intern('m')+sum(intern('ordnung')))+intern('parameter')
             u(j,1)=1;
         end
         ret=u;
     case 'ordnung'
         ret=[2  1  1];
     case 'parameter'
         ret=0;
     case 'c'
         ret=[];
     case 'n'
         ret=3;
     case 'Infinite'
         ret=false;
     case 'EVP'
         ret=true;
     case 'Endpoint'
         ret=[];
     case 'standard'
         ret=true;
     case 'rho'
         ret=[1 8];
     case 'x1'
         ret=1:0.01:2;
     case 'g'
         ret=[z(1,3)-((1+z(1,2)^2)/(z(1,1)-lamda));;z(2,2)-(0);z(3,2)-(z(1,1)^2);];
     case 'Dg'
         switch D
