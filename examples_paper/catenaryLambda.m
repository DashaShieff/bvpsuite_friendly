function [ret]=catenaryLambda(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=1;
     case 'c'
         ret=[0  1];
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
         ret=[1 5];
     case 'x1'
         ret=0:0.1:2;
     case 'g'
         ret=[0-(z(1,3)-(1+z(1,2)^2)/(z(1,1)-p(1)))];
     case 'Dg'
         switch D
          case 1
                 ret=[-(z(1,2)^2+1)/(p(1)-z(1,1))^2 ;];
          case 2
                 ret=[-(2*z(1,2))/(p(1)-z(1,1)) ;];
          case 3
                 ret=[-1 ;];
         end
     case 'Dpg'
         ret=[(z(1,2)^2+1)/(p(1)-z(1,1))^2 ;];
                 %Additional conditions
                 %zc(Komponente,Ableitung,Intervallstelle c_i)
     case 'R'
         ret=[zc(1,2,1)-(0);zc(1,1,2)-(1);p(1)-(0)];
     case 'DR'
          switch R
        case 1
            switch D
                      case 'c1_1'
                          ret=[0 ;];
                      case 'c1_2'
                          ret=[1 ;];
                      case 'c2_1'
                          ret=[0 ;];
                      case 'c2_2'
                          ret=[0 ;];
                  end
        case 2
            switch D
                      case 'c1_1'
                          ret=[0 ;];
                      case 'c1_2'
                          ret=[0 ;];
                      case 'c2_1'
                          ret=[1 ;];
                      case 'c2_2'
                          ret=[0 ;];
                  end
        case 3
            switch D
                      case 'c1_1'
                          ret=[0 ;];
                      case 'c1_2'
                          ret=[0 ;];
                      case 'c2_1'
                          ret=[0 ;];
                      case 'c2_2'
                          ret=[0 ;];
                  end
          end
     case 'DpR'
         ret=[0 ;0 ;1 ;];
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
%[2]%#dim
%parameter
%1%#parameter
%c
%[0 1]%#c
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
%0:0.1:2%#x1
%abstol
%1e-06%#abstol
%reltol
%0.001%#reltol
%auto
%false%#auto
%g
%0=z1''-(1+z1'^2)/(z1-p1)%#g
%r1
%z1'(c1)=0;z1(c2)=1;p1=0%#r1
%startprofilgitter
%#startprofilgitter
%startprofilwerte
%#startprofilwerte
