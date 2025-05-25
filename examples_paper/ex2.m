function [ret]=ex2(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=false;
     case 'EVP'
         ret=false;
     case 'Endpoint'
         ret=[];
     case 'standard'
         ret=true;
     case 'rho'
         ret=[1 2];
     case 'x1'
         ret=0:0.5:1;
     case 'g'
         ret=[z(1,3)+3/t*z(1,2)-(-(81)*z(2,1)-2*(1000)+z(1,1)*z(2,1));z(2,3)+3/t*z(2,2)-((81)*z(1,1)-(1/2)*z(1,1)^2)];
     case 'Dg'
         switch D
          case 1
                 ret=[-z(2,1) 81-z(1,1) ;z(1,1)-81 0 ;];
          case 2
                 ret=[3/t 0 ;0 3/t ;];
          case 3
                 ret=[1 0 ;0 1 ;];
         end
     case 'Dpg'
         ret=[0;0;];
                 %Additional conditions
                 %za(Komponente,Ableitung)
     case 'R'
         ret=[za(1,2)-(0);za(2,2)-(0);zb(1,1)-(0);zb(2,2)+(1-1/3)*zb(2,1)-(0)];
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
                          ret=[0 0 ;];
                      case 'a2'
                          ret=[0 1 ;];
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
                          ret=[1 0 ;];
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
                          ret=[0 2/3 ;];
                      case 'b2'
                          ret=[0 1 ;];
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
%%#Endpoint
%startwert
%1%#startwert
%dim
%[2 2]%#dim
%parameter
%0%#parameter
%c
%[]%#c
%variablen
%musquare=81;gamma=1000;%#variablen
%standard
%true%#standard
%Infinite
%false%#Infinite
%EVP
%false%#EVP
%rho
%[1 2]%#rho
%x1
%0:0.5:1%#x1
%abstol
%0.0001%#abstol
%reltol
%0.0001%#reltol
%auto
%false%#auto
%g
%z1''+3/t*z1'=-musquare*z2-2*gamma+z1*z2;  z2''+3/t*z2'=musquare*z1-(1/2)*z1^2%#g
%r1
%z1'(a)=0;z2'(a)=0;z1(b)=0;z2'(b)+(1-1/3)*z2(b)=0%#r1
%startprofilgitter
%#startprofilgitter
%startprofilwerte
%#startprofilwerte
