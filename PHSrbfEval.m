function AP = PHSrbfEval(Epoints,Ctrs,ScalePar)
% This function computes the evaluation RBF matrix at Gaussian points 

% Inputs:
% Epoints: evalution points (Gaussian points)
% Ctrs: interpolation points (center of triangles)
% ScalePar: scaling parameter 

% Outputs:
% AP: Evaluation RBF matrix 
%                      
xc = Ctrs(1,:);                    
Ctrs = (Ctrs-xc)/ScalePar;          
Epoints = (Epoints-xc)/ScalePar;
rbfpar = 2; polyorder = 2;
A = KerMat(Epoints,Ctrs,'1','tp',rbfpar);
P = PolyMat(Epoints,'1',polyorder); 
AP = [A P];