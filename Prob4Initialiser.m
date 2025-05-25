function smesh = Prob4Initialiser(xmesh,sineScale,y1,t)
    a = sineScale.*sin(pi.*t)./(2.*(xmesh(1)-xmesh(end)));
    b = -2*a*xmesh(end);
    %c = (y2 - a.*xmesh(end).^2 - b.*xmesh(end));
    c = (y1 - a.*xmesh(1).^2 - b.*xmesh(1));
    smesh = a.*xmesh.^2 + b.*xmesh + c.*ones(1,length(xmesh));
end

