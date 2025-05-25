function [smeshNew, aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto]...
    = GenerateProb2ODESysMesh(xmesh, smeshOld, calc_aw, lambda, aL) 
% Sets up the ODE system the full version of the Problem 2 ODE. This
% involes generating the 's' mesh and the 'antwort' settings for the BVP
% solver.

%%   Generate x axis mesh and s axis mesh for truncated ODE system

% Perform integration and save as 's' mesh matrix
smesh1 = smeshOld(1,:);
smeshd = gradient(smesh1,xmesh);
smeshintegral = (1 + smeshd.^2).^(1/2);
smesh2 = cumtrapz(xmesh,smeshintegral);
smesh3 = lambda.*ones(1,length(xmesh));
smeshNew = [smesh1;smesh2;smesh3];

%% Define 'antwort' settings 
[aw1, aw2, aw3, aw4, aw5, aw6, aw7, awEVP, awInf, awEP, awAuto] ...
    = deal(0,0,0,0,0,0,0,0,0,0,0);
   

%initialise numerator in the first ode:
%vel is symmetric
dX = xmesh(end) - xmesh(1);
dY = smeshNew(1,end) - smeshNew(1,1);
p1p2Dist = sqrt(dX^2 + dY^2);
Px = dX/p1p2Dist;
Py = dY/p1p2Dist;

numerator = strcat("(-(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*z1'^7" + ...
"+( -2*(2.7925)*" + num2str(Px^2) + ") *z1'^6" +...
"+(-(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + "))*z1'^5" + ...
"+(-5*(2.7925)*" + num2str(Px^2) + " - (2.7925)*" + num2str(Py^2) + ")*z1'^4" + ...
"+(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*z1'^3" + ...
"+(-4*(2.7925)*" + num2str(Px^2) + " -2*(2.7925)*" + num2str(Py^2) + ")*z1'^2" + ...
"+(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*z1'" + ...
"-(2.7925)*" + num2str(Px^2) + " - (2.7925)*" + num2str(Py^2) + ")");


if isstring(lambda) == 0
    %intialise denominator terms in the first ode:
    
    denominator = strcat("((-2*(2.7925)*" + num2str(Px^2) + "*t - 2*(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*z1 )*z1'^5" + ...
    "+(-5*(2.7925)*" + num2str(Px^2) + "*t - 5*(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*z1 )*z1'^3" + ...
    "+(-" + num2str(lambda) +" + (2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*t + (2.7925)*" + num2str(Py^2) + "*z1)*z1'^2" + ...
    "+(-3*(2.7925)*" + num2str(Px^2) + "*t - 3*(2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*z1)*z1'" + ...
    "+(-" + num2str(lambda) +" + (2.7925)*(" + num2str(Px) + ")*(" + num2str(Py) + ")*t + (2.7925)*" + num2str(Py^2) + "*z1))");

    %initialise other ode equations:
    odeone = strcat("z1'' =" + numerator + "/" + denominator);
    odetwo = "z2' = sqrt(1 + z1'^2)";
    odethree = "z3' = 0"; 
    
    %odetwo = "z1' = z2";
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
    aw6 = strcat("[" + odeone + ";"  + odetwo + ";" + odethree + "]");
end

if calc_aw(7) % antwort7, Define Boundary Condition
    %sinVal = 2*pi*timePoint;num2str(cos(pi*timePoint))
    %boundaryFree = sqrt((1.2919 - 2*2.7925 - lambda)/(2.7925 + lambda));
    aw7 = strcat("[z1(a) ="+ smeshNew(1,1) + "; z1'(b) = "+ 0 +  "; z2(a) = "+ 0 + "; z2(b) = "+ aL + "]");
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

