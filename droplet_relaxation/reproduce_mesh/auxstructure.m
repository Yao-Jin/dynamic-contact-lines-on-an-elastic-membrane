function [geom] = auxstructure(node,elem,Bx,By)

% construct auxiliary data structures for 2D triangulation
% t2v,neighbor,elem2edge,edge2elem,edge,bdEdge,node2edge,e2v
N = size(node,1);   NT = size(elem,1);

% t2v(t,i)==1 means the i-th node is a vertex of triangle t
t2v = sparse([1:NT,1:NT,1:NT],elem,1,NT,N);

totalEdge = sort([elem(:,[2,3]);elem(:,[3,1]);elem(:,[1,2])],2);
[edge,i2,j] = unique(totalEdge,'rows','legacy');
elem2edge = reshape(j,NT,3); % local index of edges in each triangle to its global index
i1(j(3*NT:-1:1)) = 3*NT:-1:1; i1 = i1';
k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
k2 = ceil(i2/NT); t2 = i2 - NT*(k2-1);  
ix = (i1~=i2);
% indices of neighbor triangles of a triangle
neighbor = accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT,3]);
% edge2elem(k,1:2) indices of two triangles sharing edge k
% edge2elem(k,3:4) local indices of the edge k in two triangles, elem2edge(edge2elem(k,1),edge2elem(k,3)) = k
edge2elem = [t1,t2,k1,k2];
% edges and their indices on the bd
bdEdge = edge((i1==i2),:);
bdEdgeind = find(i1==i2);

NE = size(edge,1);
% node2edge(i,j) = k means edge(k,:) = [i,j] or [j,i]
node2edge = sparse(edge(:,[1,2]),edge(:,[2,1]),[1:NE,1:NE],N,N);
% e2v(e, v) = 1 if v is a vertex of e
e2v = sparse([1:NE,1:NE],[edge(:,1);edge(:,2)],1,NE,N);

% judge nodes on boundary or not
bdnode = zeros(N,1);
bdnode(unique(reshape(bdEdge,2*size(bdEdge,1),1))) = 1;

geom = [];
geom.node = node;
geom.elem = elem;
geom.t2v = t2v;
geom.neighbor = neighbor;
geom.elem2edge = elem2edge;
geom.edge2elem = edge2elem;
geom.edge = edge;
geom.bdEdge = bdEdge;
geom.bdEdgeind = bdEdgeind;
geom.node2edge = node2edge;
geom.e2v = e2v;
geom.bdnode = bdnode; % judge nodes on boundary or not



% identify indices of boundary edges on each side
bdedgeup = zeros(size(edge,1),1);
bdedgeleft = zeros(size(edge,1),1);
bdedgeright = zeros(size(edge,1),1);
bdedgedown = zeros(size(edge,1),1);
for i=1:size(bdEdgeind,1)
    if (abs(node(edge(bdEdgeind(i),1),2)-By)<1e-8) && (abs(node(edge(bdEdgeind(i),2),2)-By)<1e-8)
        bdedgeup(bdEdgeind(i)) = 1;
        continue;
    end    
    if (abs(node(edge(bdEdgeind(i),1),1)+Bx)<1e-8) && (abs(node(edge(bdEdgeind(i),2),1)+Bx)<1e-8)
        bdedgeleft(bdEdgeind(i)) = 1;
        continue;
    end
    if (abs(node(edge(bdEdgeind(i),1),1)-Bx)<1e-8) && (abs(node(edge(bdEdgeind(i),2),1)-Bx)<1e-8)
        bdedgeright(bdEdgeind(i)) = 1;
        continue;
    end
    bdedgedown(bdEdgeind(i)) = 1;
end

geom.bdedgeup = bdedgeup;
geom.bdedgedown = bdedgedown;
geom.bdedgeleft = bdedgeleft;
geom.bdedgeright = bdedgeright;

% identify indices of boundary nodes on each side
bdnodeup = zeros(size(node,1),1);
bdnodedown = zeros(size(node,1),1);
bdnodeleft = zeros(size(node,1),1);
bdnoderight = zeros(size(node,1),1);
temp = edge(bdedgeup(:)==1,1:2);
bdnodeup(unique(reshape(temp,2*size(temp,1),1))) = 1;
temp = edge(bdedgedown(:)==1,1:2);
bdnodedown(unique(reshape(temp,2*size(temp,1),1))) = 1;
temp = edge(bdedgeleft(:)==1,1:2);
bdnodeleft(unique(reshape(temp,2*size(temp,1),1))) = 1;
temp = edge(bdedgeright(:)==1,1:2);
bdnoderight(unique(reshape(temp,2*size(temp,1),1))) = 1;

% exclude the overlaped nodes to upbd on the lrbd
bdnodeleft(bdnodeup(:)==1)=0;
%bdnodeleft(bdnodedown(:)==1)=0;
bdnoderight(bdnodeup(:)==1)=0;
%bdnoderight(bdnodedown(:)==1)=0;

geom.bdnodeup = bdnodeup;
geom.bdnodedown = bdnodedown;
geom.bdnodeleft = bdnodeleft;
geom.bdnoderight = bdnoderight;

end