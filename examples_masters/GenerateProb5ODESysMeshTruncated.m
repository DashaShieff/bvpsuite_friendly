function [smeshNew, aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto]...
    = GenerateProb5ODESysMeshTruncated(xmesh, smeshOld, calc_aw, lambda, aL, slope) 
% Sets up the ODE system for a truncated version of the Problem 3 ODE. This
% involes generating the 's' mesh and the 'antwort' settings for the BVP
% solver.

%%   Generate x axis mesh and s axis mesh for truncated ODE system

% Data rearragement 
smeshd = gradient(smeshOld(1,:),xmesh); % Generate s'(x)
smeshdd = gradient(smeshd,xmesh);
smeshddd = gradient(smeshdd,xmesh);
smeshintegral = ((1.2919).*smeshd - (2.7925).*smeshd)./((1 + smeshd.^2));
smeshintegral = (1 + smeshd.^2).^(1/2);
% Perform integration and save as 's' mesh matrix
smesh1 = smeshOld(1,:);
smesh2 = cumtrapz(xmesh,smeshintegral);
%smesh2 = ((1.2919) - (2.7925)).*(xmesh - xmesh(1))./2;
%smesh2 = smeshddd;
smesh3 = smeshOld(1,:);
smesh4 = smeshd;
smesh5 = smeshdd;
lambdamesh = lambda*ones(1,size(xmesh,2));
%smeshNew = [smesh1;smesh2;smesh3;smesh4;smesh5;lambdamesh];
smeshNew = [smesh1;smesh2;lambdamesh];
%smeshNew = [smesh1];
%% Define 'antwort' settings 
[aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto] ...
    = deal(0,0,0,0,0,0,0,0,0,0,0);
    
if isstring(lambda) == 0
    %intialise terms in the first ode:
    %sdthree = "2*(2.7925)*(1 + z1'^2)^2 ";
    %sdtwo = strcat("((1.2919) + "+num2str(lambda)+"*((2.7925)-3*(1.2919))*z1' + 5*(2.7925)*z1'^2)*(1 + z1'^2)");
    %sdone = strcat("3*(1.2919)*z1'^2 + 3*(2.7925)*z1'^4 + 3*"+num2str(lambda)+"*(1.2919)*z1'^3 - 3*"+num2str(lambda)+"*(2.7925)*z1'^3");
    %sdthree = strcat("("+num2str(lambda)+" + 2*(2.7925))*(1 + z1'^2)^2 ");
    %sdtwo = strcat("((1.2919) + "+num2str(lambda)+"*z1'^2 + 5*(2.7925)*z1'^2)*(1 + z1'^2)");
    %sdone = strcat("3*(1.2919)*z1'^2 + 3*(2.7925)*z1'^4 ");
    sdthree = strcat("(z3 + 2*(2.7925))*(1 + z1'^2)^2 ");
    sdtwo = strcat("((1.2919) + z3*z1'^2 + 5*(2.7925)*z1'^2)*(1 + z1'^2)");
    sdone = strcat("3*(1.2919)*z1'^2 + 3*(2.7925)*z1'^4 ");
    %sdtwo = strcat("((1.2919) + z3*((2.7925)-3*(1.2919))*z1' + 5*(2.7925)*z1'^2)*(1 + z1'^2)");
    %sdone = strcat("3*(1.2919)*z1'^2 + 3*(2.7925)*z1'^4 + 3*z3*(1.2919)*z1'^3 - 3*z3*(2.7925)*z1'^3");
    %sdtwo = strcat("((1.2919) + 5*(2.7925)*z1'^2)*(1 + z1'^2)");
    %sdone = strcat("3*(1.2919)*z1'^2 + 3*(2.7925)*z1'^4 ");
    %initialise other ode equations:
    odeone = strcat("0 = z1''*(" + sdthree + "-" + sdtwo + "+" + sdone + ")");
    %odeone = "z2' = 0";
    %odetwo = "z2' = ((1.2919)*z1' - (2.7925)*z1')/(1 + z1'^2)";
    odetwo = "z2' = (1 + z1'^2).^(1/2)";
    %odetwo = "z4' = 1";
    odethree = "z1' = z2";
    odefour = "z2' = z3";
    odefive = "z3' = z4";
    odesix = "z3' = 0";
end

if calc_aw(1) % antwort1, not sure
    aw1 = '1'; 
end

if calc_aw(2) % antwort2, Define Order of ODE System
    aw2 = '[2,1,1]';
end

if calc_aw(3) % antwort3, automatically find how many collocation points are needed
    aw3 = 'false';     
end

if calc_aw(4) % antwort4, rho based on number of collocation points (starting coefficient matrix)
    aw4 = '[1/2]'; 
end

if calc_aw(5) % antwort5, x mesh
    aw5 = strcat('[',num2str(xmesh),']'); 
end

if calc_aw(6) && (isstring(lambda) == 0) % antwort6, Define ODE System
    %aw6 = strcat("["+ odeone + ";"  + odetwo + ";" + ...
    %        odethree + ";" + odefour + ";" + odefive + ";" + odesix + "]");
    aw6 = strcat("["+ odeone + ";"  + odetwo + ";" + odesix +"]");
    %aw6 = strcat("["+ odeone +"]");
end

if calc_aw(7) % antwort7, Define Boundary Condition
    %integralx1 = ((1.2919).*smeshd - (2.7925).*smeshd)./((1 + smeshd.^2));
    integralx1 = (1 + smeshd.^2).^(1/2);
    Stress = trapz(xmesh,integralx1);
    Stress = sqrt(smeshNew(1,end)^2 + smeshNew(1,1)^2);
    Stress = ((1.2919) - (2.7925))/2;
    %sinVal = 2*pi*timePoint;z1'(b) ="+ num2str(cos(pi*slope))num2str(slope)
    %aw7 = strcat("[z1(b) ="+ smeshNew(1,end)+ ";z1(a) = "+ smeshNew(1,1)+";z3(b) = "+ 0 +";z4(b) = "+ 0 +" ]");
    %aw7 = strcat("[z2(a) ="+ num2str(slope) + ";z1(a) = "+ smeshNew(1,1)+";z3(b) = "+ 0 +";z4(b) = "+ 0 + ";z4(a) = "+ 0 +" ]");
    aw7 = strcat("[z1(b) ="+ smeshNew(1,end) + ";z1(a) = "+ smeshNew(1,1)+";z2(b) = "+ num2str(Stress) + ";z2(a) = "+ 0 + " ]");
    %aw7 = strcat("[z1(b) ="+smeshNew(1,end) + ";z1(a) = "+ smeshNew(1,1)+" ]");
end

if calc_aw(8) % antwort_EVP, is it an eigenvalue problem?
    awEVP = 'false';
end

if calc_aw(9) % antwort_Infinite, does the x mesh lead to infinity?
    awInf = 'false';
end

if calc_aw(10) % antwort_endpoint, infinity endpoint (unsure)
    awEP = [];
end

if calc_aw(11) % antwort_auto, not sure
    awAuto = 'false';
end


end

