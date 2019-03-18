function [ matrixStandard, CMean, CStd] = Standard(matrixIn, dim, CMean, CStd)
%%  Standard
%
%   Standardization of input matrix in specified dimension.
%   This standardization is suitable for rescaleing the dynamic range of
%   vector or matrix.
%   This function is able to compute Mean and Std vectors or
%   use precomputed vectors. This is suitable for standardization of Train
%   set and use the Mean and Std vectors from the Train set to perform 
%   standardization on Validation and Test sets.
%
%
%   Syntax: [ matrixStandard, CMean, CStd] = Standard(matrixIn, dim, CMean, CStd)
%
%
%       matrixIn:  
%                   Input matrix or vector.
%       
%       dim:        
%                   Desired dimention in which the standardization will be
%       performed. 
%                   dim = 1 : Standardization by row.
%                   dim = 2 : Standardization by column.
%
%       CMean:      
%                   Input/Output vector of mean values.
%
%       CStd:       
%                   Input/Output vector of std values.
%   
%       matrixStandard:
%                   
%                   Output standardized matrix or vector. 

%%  Input Mean and Std vectors.
%   In case there are no input Mean and Std vectors compute their values.
       
        switch nargin
              
            case 2
                
                CMean = mean ( matrixIn, dim );
                
                CStd  = std ( matrixIn, [], dim );
                
        end
%%  Prepare matrix.
%   
%   Prepare matrix full of ones with the needed dimensions.
    [ d1, d2 ] =   size ( matrixIn );
    
    ss         =   ones ( d1, d2);

%%  Compute the output standardized matrix.
%
%   The standardization is computed by the means of the following formula.
%
%   Xs = (X - Mean)/Std
%   
%   For each element of the matrix subtract the Mean and divide by the
%   Std elemetwise.
%

matrixStandard = (matrixIn - gmultiply(ss, CMean)) ./ gmultiply(ss, CStd);

end

