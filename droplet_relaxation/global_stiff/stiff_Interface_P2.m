function [stiff_interface] = stiff_Interface_P2(geom,para,sln)

% construct stiffness matrices for interface elements

%% ----------------------------------------------------------------------------
    % read geometry data
    N = size(geom.node,1); Ne = size(geom.edge,1);
    ffi2node = geom.ffi2node; ffi2edge = geom.ffi2edge;
    Nffi = para.Nffi; pffi_nodal = sln.pffi(1:Nffi,1:2);
    
%% ---------------------------------------------------------------------------
    % read model parameters
    Ca = para.Ca;
    
%% ---------------------------------------------------------------------------
    % initialize indices for element 
    linelemP2 = [1:Nffi-1; Nffi+1:2*Nffi-1; 2:Nffi]';
    elemP2online = [ffi2node(1:end-1)' ffi2edge(1:end)+N ffi2node(2:end)'];
    
    i1 = zeros(9*(Nffi-1),1); j1 = zeros(9*(Nffi-1),1); s1 = zeros(9*(Nffi-1),1); 
    i2 = zeros(9*(Nffi-1),1); j2 = zeros(9*(Nffi-1),1); s2 = zeros(9*(Nffi-1),1);
    i3 = zeros(9*(Nffi-1),1); j3 = zeros(9*(Nffi-1),1); s3 = zeros(9*(Nffi-1),1); 
        
    ind = 0;
    for t = 1:Nffi-1
        len = sqrt((pffi_nodal(t+1,1)-pffi_nodal(t,1))^2+(pffi_nodal(t+1,2)-pffi_nodal(t,2))^2);
        Ls1 = 1/Ca/len*getStiffonRef1DP2der;
        Ls2 = len*getStiffonRef1DP2P2;
        Ls3 = len*getStiffonRef1DP2P2;
 
        for ti = 1:3
            for tj = 1:3
                ind = ind + 1;
                i1(ind) = elemP2online(t,ti); j1(ind) = linelemP2(t,tj); s1(ind) = Ls1(ti,tj);
                i2(ind) = linelemP2(t,ti); j2(ind) = linelemP2(t,tj); s2(ind) = Ls2(ti,tj);
                i3(ind) = elemP2online(t,ti); j3(ind) = linelemP2(t,tj); s3(ind) = Ls3(ti,tj);
            end
        end
        
    end
    
    stiff_interface = [];
    stiff_interface.I1 = sparse(i1,j1,s1,N+Ne,2*Nffi-1);
    stiff_interface.I2 = sparse(i2,j2,s2,2*Nffi-1,2*Nffi-1);
    stiff_interface.I3 = sparse(i3,j3,s3,N+Ne,2*Nffi-1);
    
end