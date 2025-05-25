function [ret]=ex4(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=[1  1  0  0];
     case 'parameter'
         ret=0;
     case 'c'
         ret=[];
     case 'n'
         ret=4;
     case 'Infinite'
         ret=false;
     case 'EVP'
         ret=false;
     case 'Endpoint'
         ret=[];
     case 'standard'
         ret=true;
     case 'rho'
         ret=[1 5];
     case 'x1'
         ret=0:0.1:1;
     case 'g'
         ret=[t*z(1,2)-11*z(1,1)-18*z(2,1)+3*z(3,1)-z(4,1)-(t*exp(t)*(sin(t)+cos(t))-12*exp(t)*sin(t)-15*cos(t)+15);t*z(2,2)+12*z(1,1)+19*z(2,1)-2*z(3,1)+z(4,1)-(-t*sin(t)+13*exp(t)*sin(t)+17*cos(t)-17);z(1,1)+z(2,1)+z(3,1)-(exp(t)*sin(t)+2*cos(t)-2);2*z(1,1)+3*z(2,1)+1/5*z(4,1)-(2.2*exp(t)*sin(t)+3*cos(t)-3)];
     case 'Dg'
         switch D
          case 1
                 ret=[-11 -18 3 -1 ;12 19 -2 1 ;1 1 1 0 ;2 3 0 1/5 ;];
          case 2
                 ret=[t 0 0 0 ;0 t 0 0 ;0 0 0 0 ;0 0 0 0 ;];
         end
     case 'Dpg'
         ret=[0;0;0;0;];
                 %Additional conditions
                 %za(Komponente,Ableitung)
     case 'R'
         ret=[2*za(1,1)+3*za(2,1)-(0);za(1,1)+za(2,1)-(0)];
     case 'DR'
          switch R
             case 1
                 switch D
                      case 'a1'
                          ret=[2 3 0 0 ;];
                      case 'b1'
                          ret=[0 0 0 0 ;];
                  end
             case 2
                 switch D
                      case 'a1'
                          ret=[1 1 0 0 ;];
                      case 'b1'
                          ret=[0 0 0 0 ;];
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
%[1 1 0 0]%#dim
%parameter
%0%#parameter
%c
%[]%#c
%variablen
%%#variablen
%standard
%true%#standard
%Infinite
%false%#Infinite
%EVP
%false%#EVP
%rho
%[1 5]%#rho
%x1
%0:0.1:1%#x1
%abstol
%1e-006%#abstol
%reltol
%0.001%#reltol
%auto
%false%#auto
%g
%t*z1'-11*z1-18*z2+3*z3-z4=t*exp(t)*(sin(t)+cos(t))-12*exp(t)*sin(t)-15*cos(t)+15;  t*z2'+12*z1+19*z2-2*z3+z4=-t*sin(t)+13*exp(t)*sin(t)+17*cos(t)-17;z1+z2+z3=exp(t)*sin(t)+2*cos(t)-2;2*z1+3*z2+1/5*z4=2.2*exp(t)*sin(t)+3*cos(t)-3%#g
%r1
%2*z1(a)+3*z2(a)=0;  z1(a)+z2(a)=0%#r1
%startprofilgitter
%#startprofilgitter
%startprofilwerte
%#startprofilwerte
