function [smeshNew, aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto]...
    = GenerateProb2ODESysMeshTruncated(xmesh, smeshOld, calc_aw, lambda, aL) 
% Sets up the ODE system for a truncated version of the Problem 2 ODE. This
% involes generating the 's' mesh and the 'antwort' settings for the BVP
% solver.

%%   Generate x axis mesh and s axis mesh for truncated ODE system

% Perform integration and save as 's' mesh matrix
smesh1 = smeshOld(1,:);
smeshNew = [smesh1];

%% Define 'antwort' settings 
[aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto] ...
    = deal(0,0,0,0,0,0,0,0,0,0,0);
    
%initialise numerator in the first ode:
%vel is x 
numerator = "((1.2919)*z1'*(1 + z1'^2) - (2.7925)*z1'*(z1' + 2)*(z1' + 1))";

if isstring(lambda) == 0
    %intialise denominator terms in the first ode:
    %vel is x 
    denominator = strcat("((2.7925)*t*(2-z1'^2) + (1.2919)*t*(2*z1'^2 - 1) +" + num2str(lambda) +"*(1 + z1'^2))");
    %initialise other ode equations:
    odeone = strcat("z1'' =" + numerator + "/" + denominator);
end

if calc_aw(1) % antwort1, not sure
    aw1 = '1'; 
end

if calc_aw(2) % antwort2, Define Order of ODE System
    aw2 = '[2]';
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
    aw6 = strcat("[" + odeone + "]");
end

if calc_aw(7) % antwort7, Define Boundary Condition
    %sinVal = 2*pi*timePoint;num2str(cos(pi*timePoint))
    %boundaryFree = sqrt((1.2919 - 2*2.7925 - lambda)/(2.7925 + lambda));
    aw7 = strcat("[z1(a) ="+ smeshNew(1,1) + "; z1(b) = "+ smeshNew(1,end) + "]");
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

