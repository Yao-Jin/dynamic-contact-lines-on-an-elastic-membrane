function [geom, para, sln] = adaptmem_even_keep(geom,para,sln)

    Nmem = para.Nmem; cl2mem1 = geom.cl2mem(1); cl2mem2 = geom.cl2mem(2); Nm = para.Nm;
    pmem = sln.pmem;
     
    %-----------------------------------------------------------------------------------------
    lenl = getlength(pmem(1:cl2mem1,:));
    lenm = getlength(pmem(cl2mem1:cl2mem2,:));
    lenr = getlength(pmem(cl2mem2:Nmem,:));
    h = lenm/(Nm-1); Nl = floor(lenl/h); Nr = floor(lenr/h);
    newcl2mem1 = Nl+1; newcl2mem2 = Nl+Nm; newNmem = Nl+Nm+Nr;
    
    newpmem = zeros(newNmem,2);    
    newpmem(1:newcl2mem1,:) = interparc(newcl2mem1,pmem(1:cl2mem1,1),pmem(1:cl2mem1,2));
    newpmem(newcl2mem1:newcl2mem2,:) = interparc(newcl2mem2-newcl2mem1+1,pmem(cl2mem1:cl2mem2,1),pmem(cl2mem1:cl2mem2,2));
    newpmem(newcl2mem2:newNmem,:) = interparc(newNmem-newcl2mem2+1,pmem(cl2mem2:Nmem,1),pmem(cl2mem2:Nmem,2));

    if para.Nmem ~= newNmem
        plot(pmem(:,1),pmem(:,2),'-ok');
        hold on;
        plot(newpmem(:,1),newpmem(:,2),'-sr');
        plot(newpmem(newcl2mem1,1),newpmem(newcl2mem1,2),'-*b');
        plot(newpmem(newcl2mem2,1),newpmem(newcl2mem2,2),'-*b');
        hold off;
    end
    
    newpmem_mid = zeros(newNmem-1,2);
    newpmem_mid(1:newNmem-1,1) = 0.5*(newpmem(1:newNmem-1,1)+newpmem(2:newNmem,1));
    newpmem_mid(1:newNmem-1,2) = spline(newpmem(1:newNmem,1),newpmem(1:newNmem,2),newpmem_mid(1:newNmem-1,1));

    sln.pmem = [newpmem; newpmem_mid];
    
    geom.cl2memold = [cl2mem1 cl2mem2];
    geom.cl2mem = [newcl2mem1 newcl2mem2];
    para.Nmem = newNmem;
    
    Nmem = para.Nmem;
    pmem = zeros(2*Nmem-1,2);
    pmem(1:2:2*Nmem-1,1:2) = sln.pmem(1:Nmem,1:2);
    pmem(2:2:2*Nmem-2,1:2) = sln.pmem(Nmem+1:2*Nmem-1,1:2);
    kmem = getcurvature(pmem); kmem(1) = sln.kmem(1); kmem(end) = sln.kmem(Nmem);
    sln.kmem = [kmem(1:2:2*Nmem-1); kmem(2:2:2*Nmem-2)];     
    
end