function [sln, para] = adaptffi_even_P2(sln,para)

    Nffil = para.Nffil; pffil = sln.pffil;
    Nffir = para.Nffir; pffir = sln.pffir;
    
    oldpffil = zeros(2*Nffil-1,2);
    oldpffil(1:2:2*Nffil-1,1:2) = pffil(1:Nffil,1:2);
    oldpffil(2:2:2*Nffil-2,1:2) = pffil(Nffil+1:2*Nffil-1,1:2);    

    oldpffir = zeros(2*Nffir-1,2);
    oldpffir(1:2:2*Nffir-1,1:2) = pffir(1:Nffir,1:2);
    oldpffir(2:2:2*Nffir-2,1:2) = pffir(Nffir+1:2*Nffir-1,1:2);      
    %-----------------------------------------------------------------------------------------
    newpffil = interparc(2*Nffil-1,oldpffil(:,1),oldpffil(:,2));
    newpffir = interparc(2*Nffir-1,oldpffir(:,1),oldpffir(:,2));
    
    sln.pffil = [newpffil(1:2:2*Nffil-1,1:2); newpffil(2:2:2*Nffil-2,1:2)]; 
    sln.pffir = [newpffir(1:2:2*Nffir-1,1:2); newpffir(2:2:2*Nffir-2,1:2)];     
end
