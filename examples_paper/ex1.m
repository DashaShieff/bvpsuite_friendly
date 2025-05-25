function [ret]=ex1(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=[2];
     case 'parameter'
         ret=0;
     case 'c'
         ret=[];
     case 'n'
         ret=1;
     case 'Infinite'
         ret=false;
     case 'EVP'
         ret=false;
     case 'Endpoint'
         ret=[];
     case 'standard'
         ret=true;
     case 'rho'
         ret=[1 8];
     case 'x1'
         ret=-1:0.02:1;
     case 'g'
         ret=[(1e-4)*z(1,3)+z(1,2)-(1+(1e-4))*z(1,1)-(0)];
     case 'Dg'
         switch D
          case 1
                 ret=[-10001/10000 ;];
          case 2
                 ret=[1 ;];
          case 3
                 ret=[1/10000 ;];
         end
     case 'Dpg'
         ret=[0;];
                 %Additional conditions
                 %za(Komponente,Ableitung)
     case 'R'
         ret=[za(1,1)-(1+exp(-2));zb(1,1)-(1+exp(-2*(1+(1e-4))/(1e-4)))];
     case 'DR'
          switch R
             case 1
                 switch D
                      case 'a1'
                          ret=[1 ;];
                      case 'a2'
                          ret=[0 ;];
                      case 'b1'
                          ret=[0 ;];
                      case 'b2'
                          ret=[0 ;];
                  end
             case 2
                 switch D
                      case 'a1'
                          ret=[0 ;];
                      case 'a2'
                          ret=[0 ;];
                      case 'b1'
                          ret=[1 ;];
                      case 'b2'
                          ret=[0 ;];
                  end
          end
     case 'DpR'
         ret=[0;0;];
     case 'm'
         if (intern('standard'))
             help=intern('rho');
             ret=help(2);
         else
             ret=length(intern('rho'));
         end
     case 'linear'
         ret=1;
 end
%Values read by bvpsuite GUI:
%Endpoint
%%#Endpoint
%startwert
%1%#startwert
%dim
%[2]%#dim
%parameter
%0%#parameter
%c
%[]%#c
%variablen
%eps=1e-4%#variablen
%standard
%true%#standard
%Infinite
%false%#Infinite
%EVP
%false%#EVP
%rho
%[1 8]%#rho
%x1
%-1:0.02:1%#x1
%abstol
%1e-006%#abstol
%reltol
%0.001%#reltol
%auto
%true%#auto
%g
%eps*z1''+z1'-(1+eps)*z1=0%#g
%r1
%z1(a)=1+exp(-2);z1(b)=1+exp(-2*(1+eps)/eps)%#r1
%startprofilgitter
%#startprofilgitter
%startprofilwerte
%#startprofilwerte
