% close all;
clear;
clc;
addpath(genpath('.'));
dbstop if error

% [0.02 0.04 0.05 0.06 0.1 0.2 0.5 1 2]
load('rectangle_wetting_data_128_ind_8.mat');

node = geom.node; elem = geom.elem;
pos = [node; 0.5*(node(geom.edge(:,1),1:2)+node(geom.edge(:,2),1:2))];
trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1));
hold on;
% quiver(node(:,1),node(:,2),sln.flu_vel(1:N,1),sln.flu_vel(1:N,2),'r');
quiver(pos(:,1),pos(:,2),sln.flu_vel(:,1),sln.flu_vel(:,2),'r');
view(2);
axis equal;
hold off;

disp(['time = ',num2str((hstry.itr-1)*para.dt)]);

node = geom.node; elem = geom.elem; N = length(node); Ne = length(elem);
elem2edge = geom.elem2edge; pmem = sln.pmemold; pffi = sln.pffiold; vel = sln.flu_vel;  
Bx = para.Bx; By = para.By;

Nx = 33; Ny = 33; hy = By/Ny;
xmem = linspace(-Bx,Bx,Nx+1);
ymem = spline(pmem(:,1),pmem(:,2),xmem);
    
cnt = 0; pos = zeros(2*Nx*Ny,2); newvel = zeros(3*Nx*Ny,2);
for i = 1:Nx+1
    x = xmem(i); y = By;
    while y>ymem(i)+hy/2
        cnt = cnt + 1; pos(cnt,1) = x; pos(cnt,2) = y;
        for ind = 1:Ne
            tx = node(elem(ind,1:3),1); ty = node(elem(ind,1:3),2);
            if inpolygon(x,y,tx,ty) == 0 
                continue; 
            end
            ref = ([tx(2)-tx(1) tx(3)-tx(1); ty(2)-ty(1) ty(3)-ty(1)])^-1*([x-tx(1);y-ty(1)]);
            X = ref(1); Y = ref(2);
            weights = [vel(elem(ind,1:3),1:2); vel(elem2edge(ind,1:3)+N,1:2)];
            bases = [(1-X-Y)*(1-2*X-2*Y) X*(2*X-1) Y*(2*Y-1) 4*X*Y 4*Y*(1-X-Y) 4*X*(1-X-Y)];
            newvel(cnt,1:2) = bases*weights;
            break;
        end
        y = y-hy;
    end
end
    
% interpolate fluid velocity on the membrane
Nmem = para.Nmem;
mem2node = geom.mem2node; mem2edge = geom.mem2edge;
for i = 1:Nx+1
    x = xmem(i); y = ymem(i);
    cnt = cnt + 1; pos(cnt,1) = x; pos(cnt,2) = y;
    for j = 1:Nmem-1
        if pmem(j,1)>x || pmem(j+1,1)<x
            continue;
        end
        S = sqrt((x-pmem(j,1))^2+(y-pmem(j,2))^2)/sqrt((pmem(j+1,1)-pmem(j,1))^2+(pmem(j+1,2)-pmem(j,2))^2);
        % P2
        weights = [vel(mem2node(j:j+1),1:2); vel(mem2edge(j)+N,1:2)];
        bases = [2*(S-1/2)*(S-1) 2*S*(S-1/2) -4*S*(S-1)];
        newvel(cnt,1:2) = bases*weights;
        break;
    end
end
% for i = 1:Nx+1
%     x = xmem(i); y = ymem(i);
%     cnt = cnt + 1; pos(cnt,1) = x; pos(cnt,2) = y;
%     for j = 1:Nmem-1
%         if pmem(j,1)>x || pmem(j+1,1)<x
%             continue;
%         end
%         S = sqrt((x-pmem(j,1))^2+(y-pmem(j,2))^2)/sqrt((pmem(j+1,1)-pmem(j,1))^2+(pmem(j+1,2)-pmem(j,2))^2);
%         if (i>=2 && i<=4) || (i>=Nx-2 && i<=Nx)
%             newvel(cnt,1:2) = vel(mem2edge(j)+N,1:2);
%             break;
%         end
%         % P2
%         weights = [vel(mem2node(j:j+1),1:2); vel(mem2edge(j)+N,1:2)];
%         bases = [2*(S-1/2)*(S-1) 2*S*(S-1/2) -4*S*(S-1)];
%         newvel(cnt,1:2) = bases*weights;
%         break;
%     end
% end   


for definecolor = 1:1
    c1 = [228 26 28]/256; c2 = [55 126 184]/256; c3 = [77 175 74]/256; c4 = [152 78 163]/256;
    c5 = [255 127 0]/256; c6 = [255 217 47]/256; c7 = [166 86 40]/256; c8 = [246 112 136]/256;
end    
    
Nmem = para.Nmem; Nffi = para.Nffi;
plot(pmem(1:Nmem,1),pmem(1:Nmem,2),'-','LineWidth',3.0,'color',c3);
hold on;
plot(pffi(1:Nffi,1),pffi(1:Nffi,2),'-r','LineWidth',3.0);
plot([-Bx -Bx],[pmem(1,2) By],'-k','LineWidth',1.5);
plot([Bx Bx],[pmem(Nmem,2) By],'-k','LineWidth',1.5);
plot([-Bx Bx],[By By],'-k','LineWidth',1.5);
quiver(pos(1:cnt,1),pos(1:cnt,2),newvel(1:cnt,1),newvel(1:cnt,2),'-','LineWidth',2,'color',c2);
axis equal;
% axis off;
xlim([-1.03,1.03]);
ylim([-0.1,1]);
% axis tight;
set(gca,'FontSize',40);
hold off;