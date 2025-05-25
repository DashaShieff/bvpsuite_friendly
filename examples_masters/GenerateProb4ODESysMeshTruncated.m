function [smeshNew, aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto]...
    = GenerateProb4ODESysMeshTruncated(xmesh, smeshOld, calc_aw, lambda, aL, slope) 
% Sets up the ODE system for a truncated version of the Problem 3 ODE. This
% involes generating the 's' mesh and the 'antwort' settings for the BVP
% solver.

%%   Generate x axis mesh and s axis mesh for truncated ODE system

% Data rearragement 
xmeshFlip = fliplr(xmesh); % Flip x for 2nd BC
smesh1Flip = fliplr(smeshOld(1,:)); % Flip s(x) for 2nd BC
smeshd = gradient(smeshOld(1,:),xmesh); % Generate s'(x)
smeshdFlip = gradient(smesh1Flip,xmeshFlip); % Flip s'(x) for 2nd BC

% Define terms for integration 
smeshintegral1 = 1./((1 + smeshd.^2).^(3/2));
smeshintegral2 = smeshd./((1 + smeshd.^2).^(1/2));
smeshintegral3 = (1 + smeshd.^2).^(1/2);
smeshintegral4 = ((2.7925).*smeshdFlip.^2 + (1.2919))./(1 + smeshdFlip.^2).^(1/2);

% Perform integration and save as 's' mesh matrix
smesh1 = smeshOld(1,:); 
smesh2 = cumtrapz(xmesh,smeshintegral1);
smesh3 = cumtrapz(xmesh,smeshintegral2);
smesh4 = cumtrapz(xmesh,smeshintegral3);
smesh5 = fliplr(cumtrapz(xmeshFlip,smeshintegral4));
smeshNew = [smesh1;smesh2;smesh3;smesh4;smesh5;];

%% Define 'antwort' settings 
[aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto] ...
    = deal(0,0,0,0,0,0,0,0,0,0,0);
    
%initialise numerator in the first ode:
numerator = "(z1'*((2.7925)*z1'^2 + (1.2919))*(1 + z1'^2)^(3/2))";

if isstring(lambda) == 0
    %intialise denominator terms in the first ode:
    sdsix = "(2.7925)*z2*z1'^6";
    sdfive = "2*(2.7925)*z3*z1'^5";
    sdfour = "((1.2919) + 2*(2.7925))*z2*z1'^4";
    sdthree = "(-2*(1.2919) + 6*(2.7925))*z3*z1'^3";
    sdtwo = strcat("((2*(1.2919) - (2.7925))*z4 + (2*(1.2919) + (2.7925))*z2 +",num2str(lambda),"- z5)*z1'^2");
    sdone = "(-2*(1.2919) + 4*(2.7925))*z3*z1'";
    sdzero = strcat("(2*(2.7925) - (1.2919))*z4 + (1.2919)*z2 +",num2str(lambda),"- z5");

    %initialise other ode equations:
    odeone = strcat("z1'' =" + numerator + "/(" + sdsix + "+" + sdfive + "+" + ...
        sdfour + "+" + sdthree + "+" + sdtwo + "+" + sdone + "+" + sdzero + ")");
    odetwo = "z2' = 1/((1 + z1'^2)^(3/2))";
    odethree = "z3' = z1'/((1 + z1'^2)^(1/2))";
    odefour = "z4' = (1 + z1'^2)^(1/2)";
    odefive = "z5' = ((2.7925)*z1'^2 + (1.2919))/(1 + z1'^2)^(1/2)";
end

if calc_aw(1) % antwort1, not sure
    aw1 = '1'; 
end

if calc_aw(2) % antwort2, Define Order of ODE System
    aw2 = '[2,1,1,1,1]';
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
    aw6 = strcat("[" + odeone + ";" + odetwo + ";" + ...
            odethree + ";" + odefour + ";" + odefive + "]");
end

if calc_aw(7) % antwort7, Define Boundary Condition
    %sinVal = 2*pi*timePoint;z1'(b) ="+ num2str(cos(pi*slope))
    aw7 = strcat("[z1(b) ="+ smeshNew(1,end) + ";z1(a) = "+ smeshNew(1,1)+ ";z2(a) = 0;z3(a) = 0; z4(a) = 0;z5(b) = 0;]");
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

