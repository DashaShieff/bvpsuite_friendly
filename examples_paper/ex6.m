function [ret]=ex6(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=[2  2];
     case 'parameter'
         ret=0;
     case 'c'
         ret=[];
     case 'n'
         ret=2;
     case 'Infinite'
         ret=true;
     case 'EVP'
         ret=false;
     case 'Endpoint'
         ret=0;
     case 'standard'
         ret=true;
     case 'rho'
         ret=[1 3];
     case 'x1'
         ret=linspace(0,1,50);
     case 'g'
         ret=[z(1,3)+2/t*z(1,2)-(4*(z(1,1)+1)*z(1,1)*(z(1,1)-0.1));(-t^2*(-2*t*z(2,2)-t^2*z(2,3)))+2/(1/t)*(-t^2*z(2,2))-(4*(z(2,1)+1)*z(2,1)*(z(2,1)-0.1));];
     case 'Dg'
         switch D
          case 1
                 ret=[-4*z(1,1)*(z(1,1)-1/10)-(4*z(1,1)+4)*(z(1,1)-1/10)-(4*z(1,1)+4)*z(1,1) 0 ;0 -4*z(2,1)*(z(2,1)-1/10)-(4*z(2,1)+4)*(z(2,1)-1/10)-(4*z(2,1)+4)*z(2,1) ;];
          case 2
                 ret=[2/t 0 ;0 0 ;];
          case 3
                 ret=[1 0 ;0 t^4 ;];
         end
     case 'Dpg'
         ret=[0;0;];
                 %Additional conditions
                 %za(Komponente,Ableitung)
     case 'R'
         ret=[za(1,2)-(0);za(2,1)-(0.1);zb(1,1)-(zb(2,1));zb(1,2)-((-1^2*zb(2,2)))];
     case 'DR'
          switch R
             case 1
                 switch D
                      case 'a1'
                          ret=[0 0 ;];
                      case 'a2'
                          ret=[1 0 ;];
                      case 'b1'
                          ret=[0 0 ;];
                      case 'b2'
                          ret=[0 0 ;];
                  end
             case 2
                 switch D
                      case 'a1'
                          ret=[0 1 ;];
                      case 'a2'
                          ret=[0 0 ;];
                      case 'b1'
                          ret=[0 0 ;];
                      case 'b2'
                          ret=[0 0 ;];
                  end
             case 3
                 switch D
                      case 'a1'
                          ret=[0 0 ;];
                      case 'a2'
                          ret=[0 0 ;];
                      case 'b1'
                          ret=[1 -1 ;];
                      case 'b2'
                          ret=[0 0 ;];
                  end
             case 4
                 switch D
                      case 'a1'
                          ret=[0 0 ;];
                      case 'a2'
                          ret=[0 0 ;];
                      case 'b1'
                          ret=[0 0 ;];
                      case 'b2'
                          ret=[1 1 ;];
                  end
          end
     case 'DpR'
         ret=[0;0;0;0;];
     case 'm'
         if (intern('standard'))
             help=intern('rho');
             ret=help(2);
         else
             ret=length(intern('rho'));
         end
     case 'linear'
         ret=0;
 end
%Values read by bvpsuite GUI:
%Endpoint
%0%#Endpoint
%startwert
%1%#startwert
%dim
%[2]%#dim
%parameter
%0%#parameter
%c
%[]%#c
%variablen
%%#variablen
%standard
%true%#standard
%Infinite
%true%#Infinite
%EVP
%false%#EVP
%rho
%[1 3]%#rho
%x1
%50%#x1
%abstol
%1e-006%#abstol
%reltol
%0.001%#reltol
%auto
%false%#auto
%g
%z1''+2/t*z1'=4*(z1+1)*z1*(z1-0.1);%#g
%r1
%z1'(a)=0; z1(b)=0.1;%#r1
%startprofilgitter
%[ 0.0225    0.1000    0.1775    0.2000    0.2225    0.3000    0.3775    0.4000    0.4225    0.5000    0.5775  0.6000    0.6225    0.7000    0.7775    0.8000    0.8225    0.9000    0.9775    1.0000    1.0231    1.1111      1.2157    1.2500    1.2862    1.4286    1.6063    1.6667    1.7317    2.0000    2.3666    2.5000    2.6493  3.3333    4.4936    5.0000    5.6351   10.0000 45]#startprofilgitter
%startprofilwerte
%[-0.3042   -0.3037   -0.3024   -0.3020   -0.3014   -0.2991   -0.2962   -0.2952   -0.2942   -0.2902   -0.2857  -0.2842   -0.2828   -0.2773   -0.2713   -0.2694   -0.2675   -0.2607   -0.2535   -0.2513   -0.2490   -0.2400  -0.2286   -0.2248   -0.2207   -0.2041   -0.1825   -0.1750   -0.1669   -0.1336   -0.0899   -0.0751   -0.0592  0.0007    0.0593    0.0729    0.0834    0.0994 0.1]#startprofilwerte
