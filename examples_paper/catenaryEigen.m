function [ret]=catenaryEigen(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)

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
         ret=[z(1,3)-((z(1,2)^2+1)/(z(1,1)-z(2,1)));z(2,2)-(0);z(3,2)-(sqrt(z(1,2)^2+1));];
     case 'Dg'
         switch D
          case 1
                 ret=[(z(1,2)^2+1)/(z(1,1)-z(2,1))^2 -(z(1,2)^2+1)/(z(1,1)-z(2,1))^2 0 ;0 0 0 ;0 0 0 ;];
          case 2
                 ret=[-(2*z(1,2))/(z(1,1)-z(2,1)) 0 0 ;0 1 0 ;-z(1,2)/(z(1,2)^2+1)^(1/2) 0 1 ;];
          case 3
                 ret=[1 0 0 ;0 0 0 ;0 0 0 ;];
         end
     case 'Dpg'
         ret=[0;0;0;];
                 %Additional conditions
                 %za(Komponente,Ableitung)
     case 'R'
         ret=[za(1,1)-(1);zb(1,1)-(2);za(3,1)-(0);zb(3,1)-(2.5);];
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
%%#variablen
%standard
%true%#standard
%Infinite
%false%#Infinite
%EVP
%true%#EVP
%rho
%[1 8]%#rho
%x1
%1:0.01:2%#x1
%abstol
%1e-06%#abstol
%reltol
%0.001%#reltol
%auto
%true%#auto
%g
%z1''=(z1'^2+1)/(z1-lambda)%#g
%r1
%z1(a)=1;z1(b)=2%#r1
%startprofilgitter
%[1      1.0101      1.0202      1.0303      1.0404      1.0505      1.0606      1.0707      1.0808      1.0909       1.101      1.1111      1.1212      1.1313      1.1414      1.1515      1.1616      1.1717      1.1818      1.1919       1.202      1.2121      1.2222      1.2323      1.2424      1.2525      1.2626      1.2727      1.2828      1.2929       1.303      1.3131      1.3232      1.3333      1.3434      1.3535      1.3636      1.3737      1.3838      1.3939       1.404      1.4141      1.4242      1.4343      1.4444      1.4545      1.4646      1.4747      1.4848      1.4949      1.5051      1.5152      1.5253      1.5354      1.5455      1.5556      1.5657      1.5758      1.5859       1.596      1.6061      1.6162      1.6263      1.6364      1.6465      1.6566      1.6667      1.6768      1.6869       1.697      1.7071      1.7172      1.7273      1.7374      1.7475      1.7576      1.7677      1.7778      1.7879       1.798      1.8081      1.8182      1.8283      1.8384      1.8485      1.8586      1.8687      1.8788      1.8889       1.899      1.9091      1.9192      1.9293      1.9394      1.9495      1.9596      1.9697      1.9798      1.9899           2] #startprofilgitter
%startprofilwerte
%[1     0.96478      0.9313     0.89949     0.86926     0.84055     0.81329     0.78742     0.76286     0.73957     0.71749     0.69656     0.67674     0.65799     0.64024     0.62347     0.60764      0.5927     0.57862     0.56537     0.55291     0.54122     0.53027     0.52003     0.51048      0.5016     0.49336     0.48574     0.47874     0.47232     0.46648      0.4612     0.45647     0.45227     0.44861     0.44546     0.44283      0.4407     0.43907     0.43793      0.4373     0.43715     0.43749     0.43833     0.43966     0.44149     0.44383     0.44667     0.45002      0.4539     0.45831     0.46326     0.46876     0.47484     0.48149     0.48874      0.4966      0.5051     0.51425     0.52407      0.5346     0.54584     0.55784     0.57061     0.58419     0.59861     0.61391     0.63012     0.64728     0.66542      0.6846     0.70486     0.72625     0.74881     0.77261     0.79769     0.82411     0.85195     0.88126     0.91212     0.94459     0.97876      1.0147      1.0525      1.0923      1.1341       1.178      1.2242      1.2728      1.3238      1.3774      1.4338       1.493      1.5552      1.6205      1.6891      1.7612      1.8369      1.9165           2] #startprofilwerte
%lambda
%-1%#lambda
