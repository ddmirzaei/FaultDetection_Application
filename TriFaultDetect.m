function IdxLapDel2 = TriFaultDetect(CntTri,Y,IndXY,IndYY,uT,PUrho,h,RBFinfo)
% This function detects the fault triangles using RBF-PU method

% Inputs:
%  CntTri: barycenter of triangles  
%  Y: patch centers 
%  IndXY: indices of centers in patches 
%  IndYY: indices of patch centers in patches 
%  uT : average values at triangles 
%  PUrho: patch centers 
%  h: meshsize
%  RBFinfo: RBF information (par, type, poly oreder, ...)

% Output:
%  IdxLapDel2: indices of fault triangules 
%
LapApp = RBF_PU(CntTri,Y,IndXY,IndYY,uT,PUrho,'L',RBFinfo);
LapNorm = abs(LapApp);
CL = 1/2; CM = 1; 
alpha2 = CL/(h^1);
IdxLapAlfa2 = (LapNorm > alpha2);
LapF = LapNorm(IdxLapAlfa2); 
med = median(LapF);
delta2 = CM*med;
IdxLapDel2 = find(LapNorm > delta2);
