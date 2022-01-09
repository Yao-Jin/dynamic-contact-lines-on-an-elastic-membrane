function [geom, para, sln, reproduce_mesh] = adaptmem_even_P2(geom,para,sln)

    Nmem = para.Nmem; cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2);
    pmem = sln.pmem(1:Nmem,:);
     
    %-----------------------------------------------------------------------------------------
    len_l = getlength(pmem(1:cl2mem1,:)); hl = len_l/(cl2mem1-1);
    len_in = getlength(pmem(cl2mem1:cl2mem2,:)); hin = len_in/(cl2mem2-cl2mem1);
    len_total = getlength(pmem); ratio = hl/hin;
    Nm = cl2mem2-cl2mem1; 
    
    reproduce_mesh = 0; newcl2mem1 = cl2mem1; newcl2mem2 = cl2mem2;
    disp(['ratio: ', num2str(ratio), ' or ', num2str(1/ratio)]);
    if ratio>1.3 || 1/ratio>1.3
        Nin = floor(len_in/len_total*(Nmem-1));
        if mod(Nin,2)==1
            Nin = Nin+1;
        end
        Nout = floor((Nmem-1-Nin)/2);
        if 2*Nout+Nin+1 == Nmem
            disp('checked');
        end
        newcl2mem1 = Nout+1; newcl2mem2 = Nout+1+Nin; 
    
        if Nm ~= Nin
            reproduce_mesh = 1;
        end
    end
    
    newpmem = zeros(Nmem,2);   mid = floor((Nmem+1)/2); 
    newpmem(1:newcl2mem1,:) = interparc(newcl2mem1,pmem(1:cl2mem1,1),pmem(1:cl2mem1,2));
    newpmem(newcl2mem1:mid,:) = interparc(mid-newcl2mem1+1,pmem(cl2mem1:mid,1),pmem(cl2mem1:mid,2));
    newpmem(mid:newcl2mem2,:) = interparc(newcl2mem2-mid+1,pmem(mid:cl2mem2,1),pmem(mid:cl2mem2,2));
    newpmem(newcl2mem2:Nmem,:) = interparc(Nmem-newcl2mem2+1,pmem(cl2mem2:Nmem,1),pmem(cl2mem2:Nmem,2));
    
%     plot(newpmem(:,1),newpmem(:,2),'o-r');
%     hold on;
%     plot(sln.pmem(:,1),sln.pmem(:,2),'s--b');
%     hold off;
    
    newpmem_mid = zeros(Nmem-1,2);
    newpmem_mid(1:Nmem-1,1) = 0.5*(newpmem(1:Nmem-1,1)+newpmem(2:Nmem,1));
    newpmem_mid(1:Nmem-1,2) = spline(newpmem(1:Nmem,1),newpmem(1:Nmem,2),newpmem_mid(1:Nmem-1,1));


    sln.pmem = [newpmem; newpmem_mid];
    geom.cl2mem = [newcl2mem1 newcl2mem2];
    
    kmem = getcurvature(sln.pmem(1:Nmem,1:2)); kmem(1) = sln.kmem(1); kmem(end) = sln.kmem(Nmem);
    sln.kmem(1:Nmem) = kmem; 
    sln.kmem(Nmem+1:2*Nmem-1) = spline(newpmem(1:Nmem,1),kmem(1:Nmem),newpmem_mid(1:Nmem-1,1));
        
end