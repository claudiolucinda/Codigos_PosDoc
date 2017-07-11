 function [beta,resid,var_cov,Jstat]=ivregression(Y,X,Z,W)
        if nargin <4
            W = (Z'*Z) \ eye(size(Z,2));
        end
        beta=(X'*Z * W * Z'*X)\(X'*Z * W * Z'*Y);
        if nargout >1
            resid=Y-X * beta;
        end
        if nargout >2
           var_cov= (X'*Z * W * Z'*X) \ eye(size(X,2));
        end
        if nargout>3
            Jstat= (Z'*Y - Z'*X*beta)' * W * (Z'*Y - Z'*X*beta);
        end
 end
