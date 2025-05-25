function [smeshNew, aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto]... 
    = GenerateProb3ODESysMesh(xmesh, smeshOld, calc_aw, lambda, aL)
% Sets up the ODE system for the Problem 3 ODE. This involes generating the 
%'s' mesh and the 'antwort' settings for the BVP solver.

%vel is symmetric
dX = xmesh(end) - xmesh(1);
dY = smeshOld(1,end) - smeshOld(1,1);
p1p2Dist = sqrt(dX^2 + dY^2);
Px = dX/p1p2Dist;
Py = dY/p1p2Dist;

%%   Generate x axis mesh and s axis mesh for ODE system
xmeshFlip = fliplr(xmesh);

smesh1 = smeshOld(1,:);
smesh1Flip = fliplr(smeshOld(1,:));

smeshd = gradient(smesh1,xmesh);
smeshdFlip = gradient(smesh1Flip,xmeshFlip);

smeshintegral1 = 1./((1 + smeshd.^2).^(3/2));
smeshintegral2 = smeshd./((1 + smeshd.^2).^(1/2));
smeshintegral3 = (1 + smeshd.^2).^(1/2);
smeshintegral4 = ((2.7925).*smeshdFlip.^2 + (1.2919))./(1 + smeshdFlip.^2).^(1/2);
smeshintegral5 = (1 + smeshd.^2).^(1/2);

smeshintegral1 = (-(2.7925)*Px.*smeshdFlip + (2.7925)*Py);
smeshintegral2 = (1 + smeshd.^2).^(1/2);

% using the initial/ calculated values of smesh to calculate the integrals
smesh2 = cumtrapz(xmesh,smeshintegral1);
smesh3 = cumtrapz(xmesh,smeshintegral2);
smesh4 = cumtrapz(xmesh,smeshintegral3);
smesh5 = fliplr(cumtrapz(xmeshFlip,smeshintegral4));
smesh6 = cumtrapz(xmesh,smeshintegral5);
lambdamesh = lambda*ones(1,size(xmesh,2));

%smeshNew = [smesh1;smesh2;smesh3;smesh4;smesh5;smesh6;lambdamesh];
smeshNew = [smesh1;smesh2;smesh3;lambdamesh];
%smeshNew = [smesh1;smesh2];
%% Define 'antwort' settings 
[aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto] ...
    = deal(0,0,0,0,0,0,0,0,0,0,0);

if isstring(lambda) == 0

    %initialise numerator in the first ode:
    %numerator = "(z1'*((2.7925)*z1'^2 + (1.2919))*(1 + z1'^2)^(3/2))";
    
    %velocity normal 
    numerator = strcat("((2.7925)*"+Px+"*(1 + z1'^2)^2 + (-(2.7925)*"+Px+"*z1' + (2.7925)*"+Py+")*z1'*(1 + z1'^2))");
    
    %denominator = strcat("(" + lambda + "+ z2)");
    denominator = strcat("(z4 + z2)");
    
    %intialise denominator terms in the first ode:
    sdsix = "(2.7925)*z2*z1'^6";
    sdfive = "2*(2.7925)*z3*z1'^5";
    sdfour = "((1.2919) + 2*(2.7925))*z2*z1'^4";
    sdthree = "(-2*(1.2919) + 6*(2.7925))*z3*z1'^3";
    sdtwo = "((2*(1.2919) - (2.7925))*z4 + (2*(1.2919) + (2.7925))*z2 + z7 - z5)*z1'^2";
    sdone = "(-2*(1.2919) + 4*(2.7925))*z3*z1'";
    sdzero = "(2*(2.7925) - (1.2919))*z4 + (1.2919)*z2 + z7 - z5";

    %initialise other ode equations:
    %odeone = strcat("z1'' =" + numerator + "/(" + sdsix + "+" + sdfive + "+" + ...
    %    sdfour + "+" + sdthree + "+" + sdtwo + "+" + sdone + "+" + sdzero + ")");
    
    odetwo = "z2' = 1/((1 + z1'^2)^(3/2))";
    odethree = "z3' = z1'/((1 + z1'^2)^(1/2))";
    odefour = "z4' = (1 + z1'^2)^(1/2)";
    odefive = "z5' = ((2.7925)*z1'^2 + (1.2919))/(1 + z1'^2)^(1/2)";
    odesix = "z6' = sqrt(1 + z1'^2)";
    odeseven = "z7' = 0";   
    
    odeone = strcat("z1'' =" + numerator + "/(" + denominator + ")");
    odetwo = "z2' = (-(2.7925)*"+Px+"*z1' + (2.7925)*"+Py+")";
    odethree = "z3' = sqrt(1 + z1'^2)";
    odefour = "z4' = 0"; 
end

if calc_aw(1) % antwort1, not sure
    aw1 = '1'; 
end

if calc_aw(2) % antwort2, Define Order of ODE System
    %aw2 = '[2,1,1,1,1,1,1]';
    aw2 = '[2,1,1,1]';
    %aw2 = '[2,1]';
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
    %aw6 = strcat("[" + odeone + ";" + odetwo +  ";" + ...
    %        odethree + ";" + odefour + ";" + odefive + ";" + odesix +...
    %        ";" + odeseven + "]");
    aw6 = strcat("[" + odeone + ";" + odetwo +  ";" + ...
            odethree + ";" + odefour + "]");
    %aw6 = strcat("[" + odeone + ";" + odetwo + "]");
end

if calc_aw(7)  && (isstring(lambda) == 0) % antwort7, Define Boundary Condition
    %aw7 = strcat('[z1(a) = 1;z1(b) = 1.1;z2(a) = 0;z3(a) = 0;'...
            %+ 'z4(a) = 0;z5(b) = 0;z6(a) = 0;z6(b) =', ...
            %num2str(aL),';z7(a) =',num2str(lambda),']');
    %aw7 = strcat("[z1(a) ="+ smeshNew(1,1)+ ";z1(b) = "+ smeshNew(1,end)+ ";z2(a) = 0;z3(a) = 0;"...
    %    + "z4(a) = 0;z5(b) = 0;z6(a) = 0;z6(b) =", ...
    %    num2str(aL),";z7(a) =",num2str(lambda),"]");
    aw7 = strcat("[z1(a) ="+ smeshNew(1,1)+ ";z1(b) = "+ smeshNew(1,end)+ ";z2(b) = 0;z3(a) = 0;z3(b) =", ...
        num2str(aL)+"]");
    %aw7 = strcat("[z1(a) ="+ smeshNew(1,1)+ ";z1(b) = "+ smeshNew(1,end)+ ";z2(b) = 0]");
    Cl = "(2.7925)/z1";
    aw7 = strcat("[z1(a) ="+ smeshNew(1,1)+ ";z1'(b) = "+ "(" + Cl + "*" + aL + ")/(1 - (" + Cl + "*" + aL + ")^2)^0.5;z2(a) = 0;z3(a) = 0;z3(b) =", ...
        num2str(aL)+"]");
    aw7 = strcat("[z1(a) ="+ smeshNew(1,1)+ ";z1'(b) = "+ smeshd(end) + ";z2(a) = 0;z3(a) = 0;z3(b) =", ...
        num2str(aL)+"]");
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


