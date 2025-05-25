clear;
%close all;
%% Input Initial Conditions
%initial guess of a straight line:
x1 = 1;
x2 = 2;
y1 = 1;
y2 = 2;
arcLength = 1.1086;%1.0125;%1.4518;%42304;%1.163257639960110;%2.316964604366866;%1.445290819416838;%1.13608;%1.01927;%1.4685;-3%1.01927;-6%1.4685;-3%1.10323;%1.86613;-3%2.27177;-3%1.39368;-100%1.13608;-6
arcLength = 1.1633;
%Timesteps = 0.23;
%Timesteps =0.419569376744;
%Timesteps = 0.5;
Slopes = 0.5;
%x2SlopeArray = zeros(1,length(Slopes));
%Timesteps = [0 0.23 0.419569376744 0.5];

%ODESystemType: 0 = problem 2 truncated, 1 = problem 2 full,
%2 = problem 3 truncated, 3 = problem 3 full, 4 = problem 4
%truncated, 5 = problem 4 full 6 = problem 2 truncated
ODESystemType = 1;
initialisation = 0; % 0 = straight line, 1 = solution as initialiser, 2 = 
% parabola that meets problem 4 boundary conditions
isX2BoundFirstOrder = 0;
Lambdas = 20:-0.1:0;
%Lambdas = -8.5;
%Lambdas = 5;
Lambdas = 13.8502;
equalSpacingInteroplation = 0:0.001:1;
if initialisation == 0
    x1tau = 1:0.01:2;
    valx1tau = linspace(y1,y2,length(x1tau));
    valx1tau = -x1tau.^2 + 4.*x1tau - 3;
elseif initialisation == 1
    solFile = 'lambda5.mat';
    %solFile = 'Inflexpoint1.02.mat';
    Path = "C:\Users\shief\Documents\Dasha's Stuff\MERSTERS\Problem Set 3\lambdaSearch\";
    Path = "C:\Users\shief\Documents\Dasha's Stuff\MERSTERS\bvpsuiteREPO\bvpsuite-edits\";
    Path = 'C:\Users\shief\Documents\Dasha''s Stuff\MERSTERS\Problem Set 2 Sym Vel\FixedEndMaximising\';
    %Path = 'C:\Users\shief\Documents\Dasha''s Stuff\MERSTERS\Problem Set 4\lambdaSearchBCsdashsatx1\slope1\';
    %Path = "C:\Users\shief\Documents\Dasha's Stuff\MERSTERS\Problem Set 2 Sym Vel\11to21.1flippedmanual\";
    %Path = "C:\Users\shief\Documents\Dasha's Stuff\MERSTERS\Problem Set 3\lambdaSearchNewVel\";
    load(strcat(Path,solFile)); 
    %y2 = valx1tau(end);
    %x1tau = flipx1tau;
    %valx1tau = flipvalx1tau;
elseif initialisation == 2
    x1tau = x1:0.0001:x2;
    valx1tau = Prob4Initialiser(x1tau,1,y1,Slopes(1));
    Px = -1;
    Py = 0;
    A = 2.7925*(Py^2 + Px^2)/(Lambdas(1));
    b = -2*A*x2;
    b = 0;
    c = y1 - (sqrt(-(A*x1 + b)^2 + 1))/A;
    valx1tau = c + (sqrt(1 - (A.*x1tau +  b).^2))./A;
end
%initial guess for smesh
xmesh = x1tau;

if ODESystemType == 0
    [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateProb2ODESysMeshTruncatedSymVel(xmesh, valx1tau(1,:), ...
            [0,0,1,1,0,0,0,0,0,0,0], "BLANK", arcLength);  
elseif ODESystemType == 1 
    [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateUHessianODE(xmesh, valx1tau(1,:), ...
            [0,0,1,1,0,0,0,0,0,0,0], Lambdas(1), arcLength);  
elseif ODESystemType == 2
    [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateProb3ODESysMeshTruncated(xmesh, valx1tau(1,:), ...
            [0,0,1,1,0,0,0,0,0,0,0], "BLANK", arcLength,y1,y2);
elseif ODESystemType == 3
    [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateProb3ODESysMesh(xmesh, valx1tau(1,:), ...
        [0,0,1,1,0,0,0,0,0,0,0], Lambdas(1), arcLength);
elseif ODESystemType == 4
    [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateProb4ODESysMeshTruncated(xmesh, valx1tau(1,:), ...
        [0,0,1,1,0,0,0,0,0,0,0], "BLANK", arcLength, Slopes(1));
elseif ODESystemType == 5
    [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateProb4ODESysMesh(xmesh, valx1tau(1,:), ...
        [0,0,1,1,0,0,0,0,0,0,0], Lambdas(1), arcLength, Slopes(1));
elseif ODESystemType == 6
     [smesh,~,~,antwort3,antwort4,~,~,~,~,~,~,~]...
    = GenerateProb5ODESysMeshTruncated(xmesh, valx1tau(1,:), ...
        [0,0,1,1,0,0,0,0,0,0,0], Lambdas(1), arcLength, Slopes(1));   
end
%figure;
%plot(xmesh,smesh(1,:),'-g');

% Writing Results to File Setup
bvpfile = 'CallBvpsuiteDummyProblem'; % File Name
bvpfilem = 'CallBvpsuiteDummyProblem.m'; % File Name with Extension
pfad = "C:\Users\shief\Documents\Dasha's Stuff\MERSTERS\bvpsuiteREPO\bvpsuite-edits\examples_paper\"; % Output Directory
parameter = '0'; % parameters for pathfollowing
c = '[]'; %any boundary points inside x mesh


% Solver Options - Numerical
standard = antwort3; % Sets collocation on/off
rho = antwort4; % Sets collocation coefficient 
abstol = 1e-3; % absolute tolerance
abstolgitter = abstol; % unsure, related to absolute tolerance
reltol = 1e-3; % relative tolerance
reltolgitter = reltol; % unsure, related to relative tolerance
maxiter = 90000; % set maximum iterations
maxfunevals = 90000; % set maximum function iterations
updatejacfactor = 0.5; % update jacobian factors
lambdamin = 0.001; % smallest damping factor
switchtoffnfactor = 0.5; % switch to FF newton factor
ausgabe = 1; % verbosity of output accuracy
TRM = 0; % use trust region method
K = 1000; % unsure, likely inital mesh size
zeichnen = 0; % unsure, likely scale factor of mesh size reduction
wiederholung = 1000; % max number of iterations
feinesgitter = 1; % set fine mesh

% Solver Options - Non-Varying Functional 
default_zf_opt = optimset('Display','off','MaxIter',maxiter,'MaxFunEvals',maxfunevals); 
defaultopt = sbvpset('AbsTol',abstol,'Basis','RungeKutta','CheckJac',0,...
    'ColPts','Equidistant','Degree',4,'DegreeSelect','auto','Display',1,...
    'fVectorized',0,'IntMaxMinRatio',10,'JacVectorized',0,'MaxMeshPts',...
    10000,'OutputFcn','','OutputSel',[],'OutputTrace',1,'RelTol',reltol,'ZfOpt',default_zf_opt);

%hold on;
Index = 0;
for Lambda = Lambdas
    for Slope = Slopes
    disp(Slope);
        Index = Index + 1;
        disp(Lambda);
        if ODESystemType == 0
            [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateProb2ODESysMeshTruncatedSymVel(xmesh, smesh, ones(1,11),...
                Lambda, arcLength); 
        elseif ODESystemType == 1
            [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateUHessianODE(xmesh, smesh, ones(1,11),...
                Lambda, arcLength); 
        elseif ODESystemType == 2
            [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateProb3ODESysMeshTruncated(xmesh, smesh, ones(1,11),...
                Lambda, arcLength,y1,y2); 
        elseif ODESystemType == 3
            [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateProb3ODESysMesh(xmesh, smesh, ones(1,11),...
                Lambda, arcLength); 
        elseif ODESystemType == 4
            [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateProb4ODESysMeshTruncated(xmesh, smesh, ones(1,11),...
                Lambda, arcLength,Slope); 
        elseif ODESystemType == 5
            [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateProb4ODESysMesh(xmesh, smesh, ones(1,11),...
                Lambda, arcLength,Slope); 
        elseif ODESystemType == 6
             [smesh, antwort1, antwort2,~,~, antwort5, antwort6, antwort7,...
                antwort_EVP, antwort_Infinite, antwort_endpoint, antwort_auto]...
            = GenerateProb5ODESysMeshTruncated(xmesh, smesh, ones(1,11),...
                Lambda, arcLength,Slope);
        end
        % Write result to File
        mfileschreiben(antwort1,antwort2,antwort3,antwort4,antwort5,antwort6,...
            antwort7,antwort_EVP,antwort_Infinite,antwort_endpoint,bvpfilem,...
            pfad,parameter,c,antwort_auto);

        rehash;
        %% Run Initial file

        % Solver Options - Varying Functional 
        %Add the bvpfile to path and test feval
        addpath(pfad);
        feval(bvpfile,'x1');

        bvpopt = sbvpset;
        bvpopt.ZfOpt = optimset(defaultopt.ZfOpt,bvpopt.ZfOpt);
        bvpopt = sbvpset(defaultopt, bvpopt);
        bvpopt.Log = 0;
        bvpopt.Private.UpdateJacFactor = updatejacfactor;
        bvpopt.Private.LambdaMin = lambdamin;
        bvpopt.Private.SwitchToFFNFactor = switchtoffnfactor;
        bvpopt.TRM = TRM;

        % Initialize Mesh (variables changed from 'bvpsuite.m)
        [x1, start]=initialmesh(bvpfile,...
            xmesh,... (neuestellen)
            smesh,... (neuewerte)
            xmesh,... (neuesx1)
            str2double(parameter),... (neuesp)
            Lambda,0); % (neueslambda)

        [~, ~, ~, x1tau, valx1tau, ~, ~, ~, ~, ~, ~, ~, ~] ....
        = meshadaptation(abstolgitter, ...
                        reltolgitter,...
                        K,...
                        bvpfile,...
                        zeichnen,...
                        x1,...
                        start,...
                        bvpopt,...
                        ausgabe,...
                        wiederholung,...
                        feinesgitter);       

        %save file of output of bvpsuite function
        EVPfileName = strcat('lambda',num2str(Lambda),'.mat');
        save(EVPfileName,'valx1tau','x1tau');
        
        interpolatedVals = interparc(equalSpacingInteroplation,x1tau,valx1tau(1,:),'linear');
        x1tau = interpolatedVals(:,1)';
        valx1tau = interpolatedVals(:,2)';
        
        % Bootstrapping
        xmesh = x1tau;
        smesh = valx1tau;
        
        %extra plots
        smeshd = gradient(smesh,xmesh);
        smeshddd = gradient(gradient(gradient(smeshd,xmesh),xmesh),xmesh);
        smeshintegral = ((1.2919).*smeshd - (2.7925).*smeshd)./((1 + smeshd.^2));

        %plot(x1tau,valx1tau(1,:),'b',x1tau,smeshddd,'r',x1tau,smeshintegral,'k');
        plot(x1tau,valx1tau(1,:));
        hold on;
        
        x2Slopes = gradient(smesh(1,:),xmesh);
        x2SlopeArray(Index) = x2Slopes(1,end);
        clear x1 x1tau valx1tau start interpolatedVals x2Slope;
    end
    
end
hold on;
%plot(x2SlopeArray);
%title('Drag Minimising s(x) Curves with increasing s''(x1) unconstrained')
%xlabel('x from x1 = 1 to x2 = 2') 
%ylabel('y from y1 = 0')
%legend({'straight line','s solution'},'Location','northeast')

function ret=mfileschreiben(antwort1,antwort2,antwort3,antwort4,antwort5,antwort6,antwort7,antwort_EVP,antwort_Infinite,antwort_endpoint,bvpfile,pfad,parameter,c,antwort_auto);
dimension=num2str(length(str2num(char(antwort2))));
ordnung=str2num(char(antwort2));
parameter=str2num(char(parameter));
c=str2num(char(c));
g=char(antwort6);
r1=char(antwort7);


if strcmp(antwort_Infinite,'true')
 
  
% antwort5 = str2num(antwort5)   
 if antwort_endpoint ~= 0
    trafo_endpoint=1/antwort_endpoint;
    %help=strcat('linspace(0,',num2str(trafo_endpoint)) ;
    %antwort5 = strcat('linspace(0,1,',antwort5,')')
     help=strcat('linspace(0,',num2str(1)) ;   
 else
    help=strcat('linspace(0,',num2str(1)) ;    
     
 end 
 antwort5 = strcat(help,',',antwort5,')');
 
 end     
    
if length(pfad)> 0
    datei=fopen(strcat(pfad,bvpfile),'w');
else
    datei=fopen(bvpfile,'w');
end
help=strrep(bvpfile,'.m','');




%C: first line in bvpfile is function definition:
fprintf(datei,'function [ret]=%s(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)\n\n',help);
fprintf(datei,'%%Inputs for bvpsuite\n%%Programmed by Georg Kitzhofer, July 2004\n\n');
fprintf(datei,'%%Kontrollnummer43753976430976655145\n');
fprintf(datei,'if (nargin<15) zc=[]; if (nargin<14) DpRpar=[]; if (nargin<13) DpRgl=[]; if (nargin<12) Dppar=[]; if (nargin<11) Dpgl=[]; if (nargin<10) p=[]; if (nargin<9) D2=[]; if (nargin<8) D1=[]; if (nargin<7) zb=[]; if (nargin<6) za=[]; if (nargin<5) z=[]; if (nargin<4) t=[]; if (nargin<3) R=[];\nif (nargin<2) D=[];end;end;end;end;end;end;end;end;end;end;end;end;end;end;');
fprintf(datei,'\nret=intern(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc);\n\nfunction ret=intern(nr,D,R,t,z,za,zb,D1,D2,p,Dpgl,Dppar,DpRgl,DpRpar,zc)\n');
fprintf(datei,'switch nr\n     case ''x0''\n         for j=1:(length(intern(''x1''))-1)*(intern(''n'')*intern(''m'')+sum(intern(''ordnung'')))+intern(''parameter'')\n');

fprintf(datei,'             u(j,1)=%s;\n',char(antwort1));

fprintf(datei,'         end\n         ret=u;\n     case ''ordnung''\n         ret=[%s];\n     case ''parameter''\n         ret=%s;\n     case ''c''\n         ret=[%s];\n     case ''n''\n         ret=%s;\n',num2str(ordnung),num2str(parameter),num2str(c),char(dimension));
    
fprintf(datei,'     case ''Infinite''\n         ret=%s;\n',char(antwort_Infinite));    
fprintf(datei,'     case ''EVP''\n         ret=%s;\n',char(antwort_EVP));
if strcmp(antwort_Infinite,'true')
 

fprintf(datei,'     case ''Endpoint''\n         ret=%s;\n',char(num2str(antwort_endpoint)));

else
      
  fprintf(datei,'     case ''Endpoint''\n         ret=[];\n');  
    
end 
fprintf(datei,'     case ''standard''\n         ret=%s;\n',char(antwort3));
fprintf(datei,'     case ''rho''\n         ret=%s;\n',char(antwort4));
fprintf(datei,'     case ''x1''\n         ret=%s;\n',char(antwort5));
fprintf(datei,'     case ''g''\n');
gausgabe=prepausgabe(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
fprintf(datei,'         ret=%s;\n',gausgabe);
fprintf(datei,'     case ''Dg''\n         switch D\n');
gsymfunktion=prepdiff(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,ordnung,parameter);
linear=1;
for oi=0:max(ordnung)
    fprintf(datei,'          case %s\n',int2str(oi+1));
  %  helpy=('');
    fprintf(datei,'                 ret=[');
    for ni1=1:str2num(char(dimension))
        for nicounter=1:str2num(char(dimension))
            respect=strcat('z',int2str(nicounter),'d',int2str(oi));
            helpx=char(diff(gsymfunktion(ni1),respect));
	        helpx = regexprep(helpx, ' ', '');
            helpx=prepausgabegdiff(helpx,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
            fprintf(datei,'%s ',helpx);
            if length(strfind(helpx,'z')) || length(strfind(helpx,'p'))>0
                linear=0;
            end
        end
        fprintf(datei,';');
    end
    fprintf(datei,'];\n');
end
fprintf(datei,'         end\n');
fprintf(datei,'     case ''Dpg''\n');
%fprintf(datei,'         switch Dpgl\n');
if length(ordnung)>0
    fprintf(datei,'         ret=[');
end
for ni=1:length(ordnung)
    for pii=1:parameter
        respect=strcat('p',int2str(pii));
        help=char(diff(gsymfunktion(ni),respect));
        help = regexprep(help, ' ', '');
        help=prepausgabegdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
        fprintf(datei,'%s ',help);
        if length(strfind(help,'z')) || length(strfind(help,'p'))>0
            linear=0;
        end
    end
    if parameter>0
        fprintf(datei,';');
    else
        fprintf(datei,'0;');
    end
end
if length(ordnung)>0
    fprintf(datei,'];');
else
    fprintf(datei,'         ret=[];');
end
fprintf(datei,'\n');


fprintf(datei,'                 %%Additional conditions\n');
if length(c)==0
    fprintf(datei,'                 %%za(Komponente,Ableitung)\n');
else
    fprintf(datei,'                 %%zc(Komponente,Ableitung,Intervallstelle c_i)\n');
end
fprintf(datei,'     case ''R''\n');%         switch R\n');
if length(c)==0  
    r1ausgabe=prepausgaber(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
   else
    
    r1ausgabe=prepausgaber_c(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
end
if length(r1ausgabe)==2
    if min(r1ausgabe=='[]')
        r1ausgabe='[]';
    end
end
fprintf(datei,'         ret=%s;\n',r1ausgabe);
fprintf(datei,'     case ''DR''\n          switch R\n');
if length(c)==0
    r1symfunktion=prepdiffr(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
   
    
    
    %Jede Randbedingung hat ihre Ableitung
    
    
    for p=1:sum(ordnung)+parameter

        fprintf(datei,'             case %i\n',p);
        fprintf(datei,'                 switch D\n');

        casew=0;
        for oi=0:max(ordnung)-1
            casew=casew+1;fprintf(datei,'                      case ''a%s''\n',int2str(oi+1));
           fprintf(datei,'                          ret=[');
            %Die �u�ere Schleife bestimmt, welche Funktionen abgeleitet werden
            for nicounter=1:str2num(char(dimension))
                
                
                respect=strcat('za',int2str(nicounter),'d',int2str(oi));
                help=char(diff(r1symfunktion(p),respect));
                help = regexprep(help, ' ', '');
                help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
                fprintf(datei,'%s ',help);
                
                if length(strfind(help,'z'))>0 || length(strfind(help,'p'))>0
                    linear=0;
                end
            end
            fprintf(datei,';');
            fprintf(datei,'];\n',int2str(oi+1));
        end
        casew=0;
        for oi=0:max(ordnung)-1
            casew=casew+1;fprintf(datei,'                      case ''b%s''\n',int2str(oi+1));
            fprintf(datei,'                          ret=[');
            for nicounter=1:str2num(char(dimension))
                respect=strcat('zb',int2str(nicounter),'d',int2str(oi));
                help=char(diff(r1symfunktion(p),respect));
		        help = regexprep(help, ' ', '');
                help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
                fprintf(datei,'%s ',help);
                if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                   linear=0;
                end
            end
            fprintf(datei,';');
            fprintf(datei,'];\n',int2str(oi+1));
        end

        fprintf(datei,'                  end\n');
    end
    fprintf(datei,'          end\n');
    fprintf(datei,'     case ''DpR''\n');
    %fprintf(datei,'         switch DpRgl\n');
    if sum(ordnung)+parameter>0
        fprintf(datei,'         ret=[');
    end
    for ni=1:sum(ordnung)+parameter
        for pii=1:parameter
            respect=strcat('p',int2str(pii));
            help=char(diff(r1symfunktion(ni),respect));
            help = regexprep(help, ' ', '');
            help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
            fprintf(datei,'%s ',help);
            if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                linear=0;
            end
        end
        if parameter>0
            fprintf(datei,';');
        else
            fprintf(datei,'0;');
        end
    end
    if sum(ordnung)+parameter>0
        fprintf(datei,'];');
    else
        fprintf(datei,'         ret=[];');
    end
    fprintf(datei,'\n');
else
    r1symfunktion=prepdiffr_c(r1,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
    %Jede Randbedingung hat ihre Ableitung
    for p=1:sum(ordnung)+parameter

        fprintf(datei,'        case %i\n',p);
        fprintf(datei,'            switch D\n');
        for ci=1:length(c)
            casew=0;
            for oi=0:max(ordnung)-1
                casew=casew+1;fprintf(datei,'                      case ''c%s_%s''\n',int2str(ci),int2str(oi+1));
                fprintf(datei,'                          ret=[');
                %Die �u�ere Schleife bestimmt, welche Funktionen abgeleitet werden
                for nicounter=1:str2num(char(dimension))
                    respect=strcat('zc',num2str(ci),'_',int2str(nicounter),'d',int2str(oi));
                    help=char(diff(r1symfunktion(p),respect));
                    help = regexprep(help, ' ', '');
                    help=prepausgaberdiff_c(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
                    fprintf(datei,'%s ',help);
                    if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                        linear=0;
                    end
                end
                fprintf(datei,';');
                fprintf(datei,'];\n',int2str(oi+1));
            end
        end
        
        
        fprintf(datei,'                  end\n');
    end
    fprintf(datei,'          end\n');
    fprintf(datei,'     case ''DpR''\n');
    %fprintf(datei,'         switch DpRgl\n');
    if sum(ordnung)+parameter>0
        fprintf(datei,'         ret=[');
    end
    for ni=1:sum(ordnung)+parameter
        for pii=1:parameter
            respect=strcat('p',int2str(pii));
            help=char(diff(r1symfunktion(ni),respect));
            help = regexprep(help, ' ', '');
            if length(c)==0
                help=prepausgaberdiff(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter);
                if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                    linear=0;
                end
            else
                help=prepausgaberdiff_c(help,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,max(ordnung),parameter,c);
                if length(strfind(help,'z')) || length(strfind(help,'p'))>0
                    linear=0;
                end
            end
            fprintf(datei,'%s ',help);
        end
        if parameter>0
            fprintf(datei,';');
        else
            fprintf(datei,'0;');
        end
    end
    if sum(ordnung)+parameter>0
        fprintf(datei,'];');
    else
        fprintf(datei,'         ret=[];');
    end
    fprintf(datei,'\n');
end


fprintf(datei,'     case ''m''\n         if (intern(''standard''))\n             help=intern(''rho'');\n');
fprintf(datei,'             ret=help(2);\n         else\n             ret=length(intern(''rho''));\n');
fprintf(datei,'         end\n     case ''linear''\n         ret=');
if linear
    fprintf(datei,'1;\n end\n');
else
    fprintf(datei,'0;\n end\n');
end
fclose(datei);

%C: ende der funktion mfileschreiben
end

% helper functions for mfileschreiben
function ret=prepausgabe(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)

for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                help=strcat('z',int2str(ni));
            case 1
                help=strcat('z',int2str(ni),'''');
            case 2
                help=strcat('z',int2str(ni),'''''');
           
        end
        helpneu=strcat('z',int2str(ni),'d',int2str(oi));
       g=strrep(g,helpneu,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
       if oi<=2      
        g=strrep(g,help,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
       end
    end
end
for oi=o:-1:0
    switch oi       
        case 2
            help='z''''';
        case 1
            help='z''';
        case 0
            help='z';
    end
    helpneu=strcat('zd',int2str(oi));
    help2=strcat('xyq(1,',int2str(oi+1),')');
    g=strrep(g,helpneu,help2);
   
   
    if oi<=2
        g=strrep(g,help,help2);
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    g=strrep(g,help,strcat('yxq(',int2str(pii),')'));
end

g=strrep(g,'xyq','z');
g=strrep(g,'yxq','p');
g=strrep(g,' ','');

%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(g,'=');
gbak=g;
for j=1:length(gleichheitszeichen)
    rest=gbak(gleichheitszeichen(j):length(gbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        g=strcat(g(1:stelle-1),')',g(stelle:length(g)));
    else
        g=strrep(g,']',')]');
    end
end
g=strrep(g,'=','-(');
%Ende Gleichheitszeichen


ret=g;
end
%%%%%%%%%
function ret=prepdiff(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,ordnung,parameter)

%Erkennen der vorkommenden z und markieren durch den String xyq...
for oi=max(ordnung):-1:0
    for ni=length(ordnung):-1:1
        switch oi
            case 0
                help=strcat('z',int2str(ni));
            case 1
                help=strcat('z',int2str(ni),'''');
            case 2
                help=strcat('z',int2str(ni),'''''');
            
        end
        helpneu=strcat('z',int2str(ni),'d',int2str(oi));
        strcat('xyq(',int2str(ni),',',int2str(oi+1),')');
        g=strrep(g,helpneu,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
       
        
         if oi<=2
            g=strrep(g,help,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
        end
    end
end


for oi=max(ordnung):-1:0
    switch oi
        case 2
            help='z''''';
        case 1
            help='z''';
        case 0
            help='z';
       
     
    end
    helpneu=strcat('zd',int2str(oi));
    help2=strcat('xyq(1,',int2str(oi+1),')');
    g=strrep(g,helpneu,help2);
    
    
    if oi<=2
        g=strrep(g,help,help2);
    end
end
%Erkennen der unbekannten Parameter p und ersetzen durch den String yxq
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    g=strrep(g,help,strcat('yxq(',int2str(pii),')'));
end
g=strrep(g,'xyq','z');
g=strrep(g,'yxq','p');
g=strrep(g,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(g,'=');
gbak=g;
for j=1:length(gleichheitszeichen)
    rest=gbak(gleichheitszeichen(j):length(gbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        g=strcat(g(1:stelle-1),')',g(stelle:length(g)));
    else
        g=strrep(g,']',')]');
    end
end
g=strrep(g,'=','-(');
%Ende Gleichheitszeichen

g=eval(strcat('inline(''',g,''',''z'',''p'',''t'')'));

z=sym(0);
for oi=0:max(ordnung)
    for ni=1:length(ordnung)
        help=strcat('z',int2str(ni),'d',int2str(oi));
        z(ni,oi+1)=str2sym(help);
    end
end
p=sym(0);
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    p(pii)=str2sym(help);
end
try
    g=g(z,p,str2sym('t'));
catch
    err('bvps_errdlg31');
    g=g(z,p,str2sym('t'));
end
ret=g;
end
%%%%%%%%%
function ret=prepausgabegdiff(g,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)
for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        help=strcat('z',int2str(ni),'d',int2str(oi));
        g=strrep(g,help,strcat('xyq(',int2str(ni),',',int2str(oi+1),')'));
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    g=strrep(g,help,strcat('yxq(',int2str(pii),')'));
end
g=strrep(g,'xyq','z');
g=strrep(g,'yxq','p');
g=strrep(g,'abs(','abs1(');
ret=g;
end
%%%%%%%%%
function ret=prepausgaber(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)

for oi=o-1:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                helpa=strcat('z',int2str(ni),'(a)');
                helpb=strcat('z',int2str(ni),'(b)');
            case 1
                helpa=strcat('z',int2str(ni),'''(a)');
                helpb=strcat('z',int2str(ni),'''(b)');
            case 2
                helpa=strcat('z',int2str(ni),'''''(a)');
                helpb=strcat('z',int2str(ni),'''''(b)');     
                
                
                
        end
        helpaneu=strcat('z',int2str(ni),'d',int2str(oi),'(a)');
        helpbneu=strcat('z',int2str(ni),'d',int2str(oi),'(b)');
        r=strrep(r,helpaneu,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
        r=strrep(r,helpbneu,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
        
        if oi<=2
            r=strrep(r,helpa,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
            r=strrep(r,helpb,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
        end
    end
end
for oi=o-1:-1:0
    switch oi
        case 1
            helpa='z''(a)';
            helpb='z''(b)';
        case 0
            helpa='z(a)';
            helpb='z(b)';
        case 2 
            helpa='z''''(a)';
            helpb='z''''(b)';
            
            
    end
    helpaneu=strcat('zd',int2str(oi),'(a)');
    helpbneu=strcat('zd',int2str(oi),'(b)');
    r=strrep(r,helpaneu,strcat('xyqa(1,',int2str(oi+1),')'));
    r=strrep(r,helpbneu,strcat('xyqb(1,',int2str(oi+1),')'));
    
    if oi<=2
        r=strrep(r,helpa,strcat('xyqa(1,',int2str(oi+1),')'));
        r=strrep(r,helpb,strcat('xyqb(1,',int2str(oi+1),')'));
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end
r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(r,'=');
rbak=r;
for j=1:length(gleichheitszeichen)
    rest=rbak(gleichheitszeichen(j):length(rbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        r=strcat(r(1:stelle-1),')',r(stelle:length(r)));
    else
        r=strrep(r,']',')]');
    end
end
r=strrep(r,'=','-(');
%Ende Gleichheitszeichen
    
    
ret=r;
end
%%%%%%%%%
function ret=prepdiffr(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)

for oi=o-1:-1:0
    for ni=str2num(char(dimension)):-1:1
        switch oi
            case 0
                helpa=strcat('z',int2str(ni),'(a)');
                helpb=strcat('z',int2str(ni),'(b)');
            case 1
                helpa=strcat('z',int2str(ni),'''(a)');
                helpb=strcat('z',int2str(ni),'''(b)');
           case 2
               helpa=strcat('z',int2str(ni),'''''(a)');
               helpb=strcat('z',int2str(ni),'''''(b)');
                
        
        end
        helpaneu=strcat('z',int2str(ni),'d',int2str(oi),'(a)');
        helpbneu=strcat('z',int2str(ni),'d',int2str(oi),'(b)');
        r=strrep(r,helpaneu,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
        r=strrep(r,helpbneu,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
      
        if oi<=2
            r=strrep(r,helpa,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
            r=strrep(r,helpb,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
        end
    end
end
for oi=o-1:-1:0
    
    
    switch oi
        case 1
            helpa='z''(a)';
            helpb='z''(b)';
        case 0
            helpa='z(a)';
            helpb='z(b)';
        case 2
            helpa='z''''(a)';
            helpb='z''''(b)';
    
    end
    helpaneu=strcat('zd',int2str(oi),'(a)');
    helpbneu=strcat('zd',int2str(oi),'(b)');
    help2a=strcat('xyqa(1,',int2str(oi+1),')');
    help2b=strcat('xyqb(1,',int2str(oi+1),')');
    r=strrep(r,helpaneu,help2a);
    r=strrep(r,helpbneu,help2b);
  
    if oi<=2
        r=strrep(r,helpa,help2a);
        r=strrep(r,helpb,help2b);
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end

r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,' ','');
%Erm�gliche Gleichheitszeichen
gleichheitszeichen=strfind(r,'=');
rbak=r;
for j=1:length(gleichheitszeichen)
    rest=rbak(gleichheitszeichen(j):length(rbak)); %extrahiere alles nach dem Gleichheitszeichen
    if length(strfind(rest,';'))>0
        stelle=strfind(rest,';');
        stelle=stelle(1)+gleichheitszeichen(j)+j-2;
        r=strcat(r(1:stelle-1),')',r(stelle:length(r)));
    else
        r=strrep(r,']',')]');
    end
end
r=strrep(r,'=','-(');
%Ende Gleichheitszeichen



r=eval(strcat('inline(''',r,''',''za'',''zb'',''p'')'));

za=sym(0);
zb=sym(0);
for oi=0:o
    for ni=1:str2num(char(dimension))
        helpa=strcat('za',int2str(ni),'d',int2str(oi));
        helpb=strcat('zb',int2str(ni),'d',int2str(oi));
        za(ni,oi+1)=str2sym(helpa);
        zb(ni,oi+1)=str2sym(helpb);
    end
end
p=sym(0);
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    p(pii)=str2sym(help);
end

try
    r=r(za,zb,p);
catch
    err('bvps_errdlg32');
    r=r(za,zb,p);
end
ret=r;
end
%%%%%%%%%
function ret=prepausgaberdiff(r,antwort1,dimension,antwort3,antwort4,antwort5,antwort6,antwort7,o,parameter)
for oi=o:-1:0
    for ni=str2num(char(dimension)):-1:1
        helpa=strcat('za',int2str(ni),'d',int2str(oi));
        helpb=strcat('zb',int2str(ni),'d',int2str(oi));
        r=strrep(r,helpa,strcat('xyqa(',int2str(ni),',',int2str(oi+1),')'));
        r=strrep(r,helpb,strcat('xyqb(',int2str(ni),',',int2str(oi+1),')'));
    end
end
for pii=parameter:-1:1
    help=strcat('p',int2str(pii));
    r=strrep(r,help,strcat('yxq(',int2str(pii),')'));
end
r=strrep(r,'xyq','z');
r=strrep(r,'yxq','p');
r=strrep(r,'abs(','abs1(');
ret=r;
end
%%%%%%%%%
%%%%%%%%%