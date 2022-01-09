function [sln, para, reproduce_mesh] = adaptffi_even_P2(sln,para)

    Nffi = para.Nffi; pffi = sln.pffi;
    oldpffi = zeros(2*Nffi-1,2);
    oldpffi(1:2:2*Nffi-1,1:2) = pffi(1:Nffi,1:2);
    oldpffi(2:2:2*Nffi-2,1:2) = pffi(Nffi+1:2*Nffi-1,1:2);    
    
    len = getlength(oldpffi); hffi = len/(Nffi-1);
    hmem = getlength(sln.pmem(1:para.Nmem,1:2))/(para.Nmem-1); 
    
    reproduce_mesh = 0;
    
    if hffi/hmem<0.8 || hffi/hmem>1/0.8
       Nffi = floor(len/hmem);  
       reproduce_mesh = 1;
    end    
    
    %-----------------------------------------------------------------------------------------
    newpffi = interparc(2*Nffi-1,oldpffi(:,1),oldpffi(:,2));
        
    sln.pffi = [newpffi(1:2:2*Nffi-1,1:2); newpffi(2:2:2*Nffi-1,1:2)]; 
    
    if para.Nffi == Nffi
        reproduce_mesh = 0;
    end
    
    para.Nffi = Nffi;

    plot(oldpffi(:,1),oldpffi(:,2),'-or');
    hold on;
    plot(newpffi(:,1),newpffi(:,2),'-sb');
    hold off;

end