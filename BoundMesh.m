function [GhostCells,ExInfo,ExTR] = BoundMesh(Info,TriCnt,h,Rect,NumBC)

% This function adds ghost cells to each side of the main mesh and generate a new extended mesh 

% Inputs:
%  Info: triangulation structure for the initial mesh
%  TriCnt: barycenter of triangles
%  h: mesh size
%  Rect: rectangular domain 
%  NumBC: number of ghost cells on each side of the domain

% Outputs:
%  GhostCells: Indices of ghost cells
%
x_left = Rect.xl; x_right = Rect.xr; y_left = Rect.yl; y_right = Rect.yr;

xr = x_right-(h*NumBC);                      
IndTr = find(TriCnt(:,1) >= xr);             
Tr = Info.Elements(:,IndTr);             
IndNr = unique(Tr);
[~,IndNrb] = find(Info.Nodes(1,:) == x_right);
[~,NumNodes] = size(Info.Nodes);
IndNGl = [];
NGl = [];
count = 0;
for i= 1:length(IndNr)
    if (any(IndNr(i) == IndNrb))
        [~,ib] = ismember(round([x_left;Info.Nodes(2,IndNr(i))]',14),round(Info.Nodes',14),'rows');
        IndNGl(1,i) = ib;
    else
        count = count+1;
        IndNGl(1,i) = NumNodes+count;
        NGl(1,count) = Info.Nodes(1,IndNr(i))-(x_right-x_left);
        NGl(2,count) = Info.Nodes(2,IndNr(i));
    end
end
TGl = changem(Tr,IndNGl,IndNr);
InfoGl.Nodes = [Info.Nodes,NGl];
InfoGl.Elements = [Info.Elements,TGl];
NumNGl = count;

% ghost cells for right side of initial mesh.
xl = x_left+(h*NumBC);
IndTl = find(TriCnt(:,1) <= xl);
Tl = Info.Elements(:,IndTl);
IndNl = unique(Tl);
[~,IndNlb] = find(Info.Nodes(1,:) == x_left);
[~,NumNodes] = size(Info.Nodes);
IndNGr = [];
NGr = [];
count = 0;
for i= 1:length(IndNl)
    if (any(IndNl(i) == IndNlb))
        [~,ib] = ismember(round([x_right;Info.Nodes(2,IndNl(i))]',15),round(Info.Nodes',15),'rows');
        IndNGr(1,i) = ib;
    else
        count = count+1;
        IndNGr(1,i) = NumNodes+NumNGl+count;
        NGr(1,count) = Info.Nodes(1,IndNl(i))+(x_right-x_left);
        NGr(2,count) = Info.Nodes(2,IndNl(i));
    end
end
TGr = changem(Tl,IndNGr,IndNl);
InfoGr.Nodes = [InfoGl.Nodes,NGr];
InfoGr.Elements = [InfoGl.Elements,TGr];
TGr = InfoGr.Elements';
PGr = InfoGr.Nodes';
TRGr = triangulation(TGr,PGr);

%  ghost cells for down side of the initial mesh.
yu = y_right-(h*NumBC);
ExCnt = incenter(TRGr);
IndTu = find(ExCnt(:,2) >= yu);
Tu = InfoGr.Elements(:,IndTu);
IndNu = unique(Tu);
[~,IndNub] = find(InfoGr.Nodes(2,:) == y_right);
[~,NumNodes] = size(InfoGr.Nodes);
IndNGd = [];
NGd = [];
count = 0;
for i = 1:length(IndNu)
    if (any(IndNu(i) == IndNub))
        [~,ib] = ismember(round([InfoGr.Nodes(1,IndNu(i));y_left]',15),round(InfoGr.Nodes',15),'rows');
        IndNGd(1,i) = ib;
    else
        count = count+1;
        IndNGd(1,i) = NumNodes+count;
        
        NGd(1,count) = InfoGr.Nodes(1,IndNu(i));
        NGd(2,count) = InfoGr.Nodes(2,IndNu(i))-(y_right-y_left);
    end
end
TGd = changem(Tu,IndNGd,IndNu);
NumNGd = count;
InfoGdo.Nodes = [InfoGr.Nodes,NGd];
InfoGdo.Elements = [InfoGr.Elements,TGd];

% ghost cells for up side of the initial mesh.
yd = y_left+(h*NumBC);
IndTd = find(ExCnt(:,2) <= yd);
Td = InfoGr.Elements(:,IndTd);
IndNd = unique(Td);
[~,IndNdb] = find(InfoGr.Nodes(2,:) == y_left);
[~,NumNodes] = size(InfoGr.Nodes);
IndNGup = [];
NGup = [];
count = 0;
for i = 1:length(IndNd)
    if (any(IndNd(i) == IndNdb))
        [~,ib] = ismember(round([InfoGr.Nodes(1,IndNd(i));y_right]',15),round(InfoGr.Nodes',15),'rows');
        IndNGup(1,i) = ib;
    else
        count = count+1;
        IndNGup(1,i) = NumNodes+NumNGd+count;
        
        NGup(1,count) = InfoGr.Nodes(1,IndNd(i));
        NGup(2,count) = InfoGr.Nodes(2,IndNd(i))+(y_right-y_left);
    end
end
TGup = changem(Td,IndNGup,IndNd);
InfoGup.Nodes = [InfoGdo.Nodes,NGup];
InfoGup.Elements = [InfoGdo.Elements,TGup];
Tnew = InfoGup.Elements';
Pnew = InfoGup.Nodes';
ExTR = triangulation(Tnew,Pnew);
GhostCells = {IndTr,IndTl,IndTu,IndTd};
ExInfo.Nodes = Pnew';
ExInfo.Elements = Tnew';
