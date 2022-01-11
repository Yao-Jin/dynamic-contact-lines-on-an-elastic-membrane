function [stiff_interface] = stiff_Interface_P2(geom,para,sln)

% construct stiffness matrices for interface elements

%% ----------------------------------------------------------------------------
    % read geometry data
    N = size(geom.node,1); Ne = size(geom.edge,1);
    ffil2node = geom.ffil2node; ffir2node = geom.ffir2node; 
    ffil2edge = geom.ffil2edge; ffir2edge = geom.ffir2edge;
    Nffil = para.Nffil; Nffir = para.Nffir; 
    pffil = sln.pffil(1:Nffil,1:2); pffir = sln.pffir(1:Nffir,1:2);
    
%% ---------------------------------------------------------------------------
    % read model parameters
    Ca = para.Ca;
    
%% ---------------------------------------------------------------------------
    % initialize indices for element 
    linelemP2 = [1:Nffil-1; Nffil+1:2*Nffil-1; 2:Nffil]';
    elemP2online = [ffil2node(1:end-1)' ffil2edge(1:end)+N ffil2node(2:end)'];
    
    i1 = zeros(9*(Nffil-1),1); j1 = zeros(9*(Nffil-1),1); s1 = zeros(9*(Nffil-1),1); 
    i2 = zeros(9*(Nffil-1),1); j2 = zeros(9*(Nffil-1),1); s2 = zeros(9*(Nffil-1),1);
    i3 = zeros(9*(Nffil-1),1); j3 = zeros(9*(Nffil-1),1); s3 = zeros(9*(Nffil-1),1); 
  
    ind = 0;
    for t = 1:Nffil-1
        len = sqrt((pffil(t+1,1)-pffil(t,1))^2+(pffil(t+1,2)-pffil(t,2))^2);
        Ls1 = 1/Ca/len*getStiffonRef1DP2der;
        Ls2 = len*getStiffonRef1DP2P2;
        Ls3 = len*getStiffonRef1DP2P2;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i1(ind) = elemP2online(t,ti); j1(ind) = linelemP2(t,tj); s1(ind) = Ls1(ti,tj);
                i2(ind) = linelemP2(t,ti); j2(ind) = linelemP2(t,tj); s2(ind) = Ls2(ti,tj);
                i3(ind) = linelemP2(t,tj); j3(ind) = elemP2online(t,ti); s3(ind) = Ls3(ti,tj);
            end
        end
    end
    
    stiff_interface = [];
    stiff_interface.Il1 = sparse(i1,j1,s1,N+Ne,2*Nffil-1);
    stiff_interface.Il2 = sparse(i2,j2,s2,2*Nffil-1,2*Nffil-1);
    stiff_interface.Il3 = sparse(i3,j3,s3,2*Nffil-1,N+Ne);
    
%% ---------------------------------------------------------------------------
    % initialize indices for element 
    linelemP2 = [1:Nffir-1; Nffir+1:2*Nffir-1; 2:Nffir]';
    elemP2online = [ffir2node(1:end-1)' ffir2edge(1:end)+N ffir2node(2:end)'];
    
    i1 = zeros(9*(Nffir-1),1); j1 = zeros(9*(Nffir-1),1); s1 = zeros(9*(Nffir-1),1); 
    i2 = zeros(9*(Nffir-1),1); j2 = zeros(9*(Nffir-1),1); s2 = zeros(9*(Nffir-1),1);
    i3 = zeros(9*(Nffir-1),1); j3 = zeros(9*(Nffir-1),1); s3 = zeros(9*(Nffir-1),1); 
  
    ind = 0;
    for t = 1:Nffir-1
        len = sqrt((pffir(t+1,1)-pffir(t,1))^2+(pffir(t+1,2)-pffir(t,2))^2);
        Ls1 = 1/Ca/len*getStiffonRef1DP2der;
        Ls2 = len*getStiffonRef1DP2P2;
        Ls3 = len*getStiffonRef1DP2P2;
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i1(ind) = elemP2online(t,ti); j1(ind) = linelemP2(t,tj); s1(ind) = Ls1(ti,tj);
                i2(ind) = linelemP2(t,ti); j2(ind) = linelemP2(t,tj); s2(ind) = Ls2(ti,tj);
                i3(ind) = linelemP2(t,tj); j3(ind) = elemP2online(t,ti); s3(ind) = Ls3(ti,tj);
            end
        end
    end
    
    stiff_interface.Ir1 = sparse(i1,j1,s1,N+Ne,2*Nffir-1);
    stiff_interface.Ir2 = sparse(i2,j2,s2,2*Nffir-1,2*Nffir-1);
    stiff_interface.Ir3 = sparse(i3,j3,s3,2*Nffir-1,N+Ne);    
end