function [TR,Info] = MeshGen(h,Rect)
% This function generates a triangular mesh on domain [xl,xr]x [yl,yr] with mesh size h

% Inputs:
% h: meshsize
% Rect: rectngular domain 

% Outputs:
% TR: a triangulation class 
% Info: information about triangulation
%
gd = [3;4;Rect.xl;Rect.xr;Rect.xr;Rect.xl;Rect.yl;Rect.yl;Rect.yr;Rect.yr];  
ns = char('rect');
ns = ns';
sf = 'rect';
dl = decsg(gd,sf,ns);
model = createpde;
geometryFromEdges(model,dl);
Info = generateMesh(model,'Hmax',h,'GeometricOrder','linear');
P = Info.Nodes; P = P';
T = Info.Elements; T = T';
TR = triangulation(T,P);
 