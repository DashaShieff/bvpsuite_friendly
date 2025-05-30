function [ret]=Problem2Explicit(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=false;
     case 'rho'
         ret=[1/2];
     case 'x1'
         ret=[1      1.0112      1.0223      1.0335      1.0446      1.0557      1.0669       1.078      1.0891      1.1002      1.1113      1.1224      1.1335      1.1445      1.1556      1.1667      1.1777      1.1887      1.1998      1.2108      1.2218      1.2328      1.2437      1.2547      1.2657      1.2766      1.2875      1.2984      1.3093      1.3202      1.3311      1.3419      1.3527      1.3635      1.3743      1.3851      1.3958      1.4065      1.4172      1.4279      1.4386      1.4492      1.4598      1.4704      1.4809      1.4914      1.5019      1.5124      1.5228      1.5332      1.5435      1.5538      1.5641      1.5743      1.5845      1.5947      1.6048      1.6149      1.6249      1.6349      1.6449      1.6548      1.6647      1.6745      1.6842      1.6939      1.7036      1.7132      1.7228      1.7323      1.7418      1.7512      1.7606      1.7699      1.7791      1.7883      1.7975      1.8065      1.8156      1.8246      1.8335      1.8424      1.8512      1.8599      1.8686      1.8773      1.8858      1.8944      1.9028      1.9112      1.9196      1.9279      1.9361      1.9443      1.9525      1.9605      1.9685      1.9765      1.9844      1.9922           2];
     case 'g'
         ret=[z(1,3)*((2.7925)*t*(2-z(1,2)^2)+(1.2919)*t*(2*z(1,2)^2-1)+z(2,1)*(1+z(1,2)^2))-((1.2919)*z(1,2)*(1+z(1,2)^2)-(2.7925)*z(1,2)*(z(1,2)^2+2)*(z(1,2)^2+1));z(2,2)-(0);z(3,2)-(sqrt(z(1,2)^2+1));];
     case 'Dg'
         switch D
          case 1
                 ret=[0 z(1,3)*(z(1,2)^2+1) 0 ;0 0 0 ;0 0 0 ;];
          case 2
                 ret=[(1117*z(1,2)^2*(z(1,2)^2+1))/200+(1117*z(1,2)^2*(z(1,2)^2+2))/200+(1117*(z(1,2)^2+1)*(z(1,2)^2+2))/400-(38757*z(1,2)^2)/10000-z(1,3)*((2087*t*z(1,2))/5000-2*z(1,2)*z(2,1))-12919/10000 0 0 ;0 1 0 ;-z(1,2)/(z(1,2)^2+1)^(1/2) 0 1 ;];
          case 3
                 ret=[z(2,1)*(z(1,2)^2+1)-(1117*t*(z(1,2)^2-2))/400+(12919*t*(2*z(1,2)^2-1))/10000 0 0 ;0 0 0 ;0 0 0 ;];
         end
     case 'Dpg'
         ret=[0;0;0;];
                 %Additional conditions
                 %za(Komponente,Ableitung)
     case 'R'
         ret=[za(1,1)-(1);zb(1,1)-(1.5);za(3,1)-(0);zb(3,1)-(1.13898);];
     case 'DR'
          switch R
             case 1
                 switch D
                      case 'a1'
                          ret=[1 0 0 ;];
                      case 'a2'
                          ret=[0 0 0 ;];
                      case 'b1'
                          ret=[0 0 0 ;];
                      case 'b2'
                          ret=[0 0 0 ;];
                  end
             case 2
                 switch D
                      case 'a1'
                          ret=[0 0 0 ;];
                      case 'a2'
                          ret=[0 0 0 ;];
                      case 'b1'
                          ret=[1 0 0 ;];
                      case 'b2'
                          ret=[0 0 0 ;];
                  end
             case 3
                 switch D
                      case 'a1'
                          ret=[0 0 1 ;];
                      case 'a2'
                          ret=[0 0 0 ;];
                      case 'b1'
                          ret=[0 0 0 ;];
                      case 'b2'
                          ret=[0 0 0 ;];
                  end
             case 4
                 switch D
                      case 'a1'
                          ret=[0 0 0 ;];
                      case 'a2'
                          ret=[0 0 0 ;];
                      case 'b1'
                          ret=[0 0 1 ;];
                      case 'b2'
                          ret=[0 0 0 ;];
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
%[2]%#dim
%parameter
%0%#parameter
%c
%[]%#c
%variablen
%Cn = 2.7925;Ct = 1.2919;%#variablen
%standard
%false%#standard
%Infinite
%false%#Infinite
%EVP
%true%#EVP
%rho
%[1/2]%#rho
%x1
%[1      1.0112      1.0223      1.0335      1.0446      1.0557      1.0669       1.078      1.0891      1.1002      1.1113      1.1224      1.1335      1.1445      1.1556      1.1667      1.1777      1.1887      1.1998      1.2108      1.2218      1.2328      1.2437      1.2547      1.2657      1.2766      1.2875      1.2984      1.3093      1.3202      1.3311      1.3419      1.3527      1.3635      1.3743      1.3851      1.3958      1.4065      1.4172      1.4279      1.4386      1.4492      1.4598      1.4704      1.4809      1.4914      1.5019      1.5124      1.5228      1.5332      1.5435      1.5538      1.5641      1.5743      1.5845      1.5947      1.6048      1.6149      1.6249      1.6349      1.6449      1.6548      1.6647      1.6745      1.6842      1.6939      1.7036      1.7132      1.7228      1.7323      1.7418      1.7512      1.7606      1.7699      1.7791      1.7883      1.7975      1.8065      1.8156      1.8246      1.8335      1.8424      1.8512      1.8599      1.8686      1.8773      1.8858      1.8944      1.9028      1.9112      1.9196      1.9279      1.9361      1.9443      1.9525      1.9605      1.9685      1.9765      1.9844      1.9922           2]%#x1
%abstol
%1e-06%#abstol
%reltol
%0.001%#reltol
%auto
%false%#auto
%g
%z1''*(Cn*t*(2-z1'^2) + Ct*t*(2*z1'^2 - 1) + lambda*(1 + z1'^2)) = Ct*z1'*(1 + z1'^2) - Cn*z1'*(z1'^2 + 2)*(z1'^2 + 1)%#g
%r1
%z1(a) = 1;z1(b) = 1.5;%#r1
%startprofilgitter
%[1      1.0112      1.0223      1.0335      1.0446      1.0557      1.0669       1.078      1.0891      1.1002      1.1113      1.1224      1.1335      1.1445      1.1556      1.1667      1.1777      1.1887      1.1998      1.2108      1.2218      1.2328      1.2437      1.2547      1.2657      1.2766      1.2875      1.2984      1.3093      1.3202      1.3311      1.3419      1.3527      1.3635      1.3743      1.3851      1.3958      1.4065      1.4172      1.4279      1.4386      1.4492      1.4598      1.4704      1.4809      1.4914      1.5019      1.5124      1.5228      1.5332      1.5435      1.5538      1.5641      1.5743      1.5845      1.5947      1.6048      1.6149      1.6249      1.6349      1.6449      1.6548      1.6647      1.6745      1.6842      1.6939      1.7036      1.7132      1.7228      1.7323      1.7418      1.7512      1.7606      1.7699      1.7791      1.7883      1.7975      1.8065      1.8156      1.8246      1.8335      1.8424      1.8512      1.8599      1.8686      1.8773      1.8858      1.8944      1.9028      1.9112      1.9196      1.9279      1.9361      1.9443      1.9525      1.9605      1.9685      1.9765      1.9844      1.9922           2]#startprofilgitter
%startprofilwerte
%[1        1.01      1.0199      1.0298      1.0395      1.0492      1.0587      1.0682      1.0776      1.0869      1.0961      1.1052      1.1142      1.1231       1.132      1.1407      1.1493      1.1578      1.1662      1.1746      1.1828      1.1909      1.1989      1.2068      1.2146      1.2223      1.2298      1.2373      1.2447      1.2519      1.2591      1.2661       1.273      1.2798      1.2865      1.2931      1.2995      1.3059      1.3121      1.3182      1.3242      1.3301      1.3358      1.3415       1.347      1.3524      1.3577      1.3629       1.368       1.373      1.3778      1.3825      1.3872      1.3917      1.3961      1.4004      1.4046      1.4087      1.4126      1.4165      1.4203      1.4239      1.4275       1.431      1.4343      1.4376      1.4408      1.4438      1.4468      1.4497      1.4525      1.4552      1.4578      1.4604      1.4628      1.4652      1.4675      1.4697      1.4718      1.4738      1.4758      1.4777      1.4795      1.4812      1.4828      1.4844      1.4859      1.4874      1.4887        1.49      1.4913      1.4924      1.4935      1.4946      1.4955      1.4964      1.4973       1.498      1.4987      1.4994         1.5]#startprofilwerte
%lambda
%-0.8%#lambda
