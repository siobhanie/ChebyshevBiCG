function [x,flag,relres,iter,resvec]=agmg(A,b,icg,tol,maxit,verbose,x0,ijob)
%   AGMG: solve the linear system of equations  'A * X = B' by means of an
%         Aggregation-based algebraic Multigrid iterative Method
%         accelerated either by the Conjugate Gradient method (CG) or by
%         the Generalized Conjugate Residual method (GCR, a variant of GMRES).
%
%   The N-by-N coefficient matrix A must be square and and sparse,
%   and the right hand side column vector B must have length N.
%   The matrix A may be real or complex; for complex matrices, however,
%   AGMG is tentative only. If A is real, B has to be real too.
%
%   X = AGMG(A,B,RESTART)
%       If RESTART==1, use CG and performs some simplifications based on
%             the assumption that the coefficient matrix A is symmetric;
%             should be used only when A is symmetric and positive definite.
%       If RESTART>=2, use GCR restarted each RESTART iterations.
%       If RESTART==0 or [], then AGMG uses the default, which is
%                     GCR restarted each 10 iterations.
%
%   X = AGMG(A,B,RESTART,TOL) specifies the tolerance of the method.
%       If TOL is [], then AGMG uses the default, 1e-6.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT) specifies the maximum number of iterations.
%       If MAXIT is [], then AGMG uses the default, 100.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE)
%       If VERBOSE==1, information is displayed on the solution process
%       If VERBOSE==0 or [], only error messages are displayed (default).
%       If VERBOSE==-1, warning messages about insufficient convergence
%                       (in the allowed maximum number of iterations)
%                       are further suppressed; i.e., AGMG works silently.
%       See the user guide provided with AGMG for a description of the
%       verbose output.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE,X0) specifies an initial guess.
%       If X0 is [], then AGMG uses the default, an all zero vector.
%
%   X = AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB)
%       If IJOB==1, perform the setup only (preprocessing:
%                   prepare all parameters for subsequent solves).
%          Then, only A, RESTART and VERBOSE input parameters are significant
%                and the calling may be: AGMG(A,[],RESTART,[],[],[],[],1).
%          The returned X is empty and other output parameters are meaningless.
%       If IJOB==2, solve only, based on previous setup.
%          Then, A may differ from the matrix supplied for set up (former
%                call with IJOB==1), but it means using a preconditioner
%                computed for a matrix to solve a system with another matrix,
%                which is not recommended in general.
%       If IJOB==202, solves only, based on previous setup.
%          The system matrix is assumed unchanged and argument A
%          is not significant.
%       If IJOB==3, the returned X is not the solution of the linear system,
%                   but the result of the action of the multigrid
%                   preconditioner on the right hand side B.
%          Then, A, TOL, MAXIT,  RESTART, and X0 are not significant;
%                the calling may be: AGMG([],B,[],[],[],[],[],3)
%          Further output parameters (besides X) are meaningless.
%       If IJOB==-1, erases the setup and releases internal memory.
%          Other input parameters are not significant and the calling may be
%                AGMG([],[],[],[],[],[],[],-1).
%          The returned X is empty and other output parameters are meaningless.
%       DEFAULT: If IJOB==0 or [],  solve the linear system (performing setup
%          and solve), and release memory.
%          (other input & ouput parameters have their usual meaning).
%
%       IJOB == 100,101,102: same as, respectively, IJOB==0,1,2,
%       but use the TRANSPOSE of the input matrix in A.
%       Hence, AGMG(A,B,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB+100)
%       is equivalent to AGMG(A',B,RESTART,TOL,MAXIT,VERBOSE,X0,IJOB),
%       but less memory consuming and often faster.
%
%   !!! IJOB==2,3,102,202 require that one has previously called AGMG
%       with IJOB==1 or IJOB==101
%
%   [X,FLAG] = AGMG(A,B,...) also returns a convergence FLAG:
%        0 if AGMG converged to the desired tolerance;
%        1 if AGMG did not reached the desired tolerance.
%
%   [X,FLAG,RELRES] = AGMG(A,B,...) also returns the relative residual
%    norm (NORM(B-A*X)/NORM(B)). If FLAG is 0, then RELRES <= TOL.
%    It may happen that RELRES > TOL while FLAG is zero when TOL is
%    below the accuracy that can be attained with a backward stable solver;
%    then RELRES is comparable to what can be obtained with a direct solver.
%    (RELRES cannot be computed when IJOB==2 since the first argument may
%     be different from the system matrix).
%
%   [X,FLAG,RELRES,ITER] = AGMG(A,B,...) also returns the number of iterations
%    actually performed.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = AGMG(A,B,...) also returns a vector of the
%    residual norms at each iteration including NORM(B-A*X0).
%
%   Example: % generate sample matrix and right hand side vector
%            % (5-point discrete Poisson operator on a 50x50 grid)
%            m=50; A=delsq(numgrid('S',m+2)); n=m*m; b=ones(n,1);
%            %
%            %solve with agmg:
%            x=AGMG(A,b,1);         % A is symmetric positive definite
%            x=AGMG(A,b);           % Consider A as a general matrix
%            %
%            x=AGMG(A,b,1,[],[],1); % Verbose output
%            % Example of tolerance below attainable accuracy:
%            [x,flag,relres,iter,resvec]=agmg(A,b,1,1e-20);
%            % AGMG report normal convergence:
%            disp('Convergence flag = '),disp(flag)
%            % The relative residual is nevertheless larger than TOL
%            % (and than reported by the verbose output):
%            disp('Relative residual = '),disp(relres)
%            % But the attained accuracy is similar to the one obtained with "\":
%            y=A\b;
%            disp('Relative residual with "\" = '),disp(norm(A*y-b)/norm(b))
%
% WARNING: this version of the AGMG function is for Matlab only;
%          Octave users: if you see this help within an Octave session,
%                        please check the location of your agmg.oct file
%
% COPYRIGHT (c) 2012-2018 Yvan Notay - ULB
% This function is part of AGMG software package, release 3.3.5
% ALL USAGE OF AGMG IS SUBJECT TO LICENSE
% Enter agmg() for detailed conditions of use.
% See the Web page <http://agmg.eu> for
%    release information, a detailed userguide and possible upgrades.
%
persistent preprocessed n notcpl
if isempty(preprocessed)
   preprocessed = 0;
end
   % Check matrix and right hand side vector inputs have appropriate size
   if (nargin < 2)
fprintf(' COPYRIGHT (c) 2012-2018 Yvan Notay - ULB \n')
fprintf(' The function agmg is part of AGMG software package \n')
fprintf(' \n')
fprintf(' ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE". \n')
fprintf(' IF YOU OBTAINED A COPY OF THIS SOFTWARE WITHOUT THIS FILE, \n')
fprintf(' PLEASE CONTACT info@agmg.eu \n')
fprintf(' \n')
fprintf(' In particular, if you have a free academic license: \n')
fprintf(' \n')
fprintf(' (1) You must be a member of an educational, academic or research institution. \n')
fprintf('     The license agreement automatically terminates once you no longer fulfill \n')
fprintf('     this requirement. \n')
fprintf(' \n')
fprintf(' (2) You are obliged to cite AGMG in any publication or report as: \n')
fprintf('     "Yvan Notay, AGMG software and documentation; \n')
fprintf('      see http://agmg.eu". \n')
fprintf(' \n')
fprintf(' (3) You may not make available to others the software in any form,  \n')
fprintf('     either as source or as a precompiled object. \n')
fprintf(' \n')
fprintf(' (4) You may not use AGMG for the benefit of any third party or for any \n')
fprintf('     commercial purposes. Note that this excludes the use within the \n')
fprintf('     framework of a contract with an industrial partner. \n')
fprintf(' \n')
fprintf(' See the Web pages <http://agmg.eu> for \n')
fprintf('    release information, a detailed userguide and possible upgrades. \n')
fprintf(' \n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
fprintf(' DICLAIMER: \n')
fprintf('    AGMG is provided on an "AS IS" basis, without any explicit or implied \n')
fprintf('    WARRANTY; see the see the file "LICENSE" for more details. \n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
fprintf('   If you use AGMG for research, please observe that your work benefits \n')
fprintf('   our past research efforts that allowed the development of AGMG. \n')
fprintf('   Hence, even if you do not see it directly, the results obtained thanks \n')
fprintf('   to the use of AGMG depend on the results in publications [1-3] below,\n')
fprintf('   where the main algorithms used in AGMG are presented and justified.   \n')
fprintf('   It is then a normal duty to cite these publications (besides citing \n')
fprintf('   AGMG itself) in any scientific work depending on the usage of AGMG, \n')
fprintf('   as you would do with any former research result you are using. \n')
fprintf(' \n')
fprintf(' [1] Y. Notay, An aggregation-based algebraic multigrid method, \n')
fprintf('    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010 \n')
fprintf('\n')
fprintf(' [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed \n')
fprintf('    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012. \n')
fprintf(' \n')
fprintf(' [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion \n')
fprintf('    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012. \n')
fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
fprintf(' \n')
fprintf(' Error using agmg \n')
fprintf(' agmg requires at least two input arguments \n')
fprintf(' enter "help agmg" for more information \n')
   return
   end
   if (nargin < 3) || isempty(icg)
      icg = 0;
   end
   if (nargin < 4) || isempty(tol)
      tol = 1e-6;
   end
   if (nargin < 5) || isempty(maxit)
      maxit = 100;
   end
   if (nargin < 6) || isempty(verbose)
      verbose=0;
   end
   if (verbose>0) verbose=6; end
   if (nargin < 7)
      x0=[];
   end
   if (nargin < 8 || isempty(ijob))
      ijob=0;
   end
   if (ijob >= 100 && ijob ~= 202)
	   ijb=ijob-100;
   else
	   ijb=ijob;
   end
   if ( (ijb < -1 || ijb > 2 || round(ijb)~=ijb) && ijob ~= 202 && ijob~=3)
          error('MATLAB:agmg:nonvalidijob',...
                'IJOB should be equal to -1, 0, 1, 2, 3, 100, 101, 102 or 202');
   elseif (ijb > 1 && ~preprocessed)
          error('MATLAB:agmg:nonvalidijob',...
                'Setup not done: IJOB should be equal to 0, 1, 100 or 101')
   elseif (ijb == -1 && ~preprocessed)
          disp('Warning: setup not done: nothing to do for IJOB == -1')
          iter=0;
          relres=[];
          x=[];
          resvec=[];
          return
   end
   if (ijb==0 || ijb==1 || ijb==2)
      if (isempty(A) ||  any(class(A)~='double') || ~issparse(A))
         error('MATLAB:agmg:NonSparseMatrix', ...
               'the input matrix A must be a sparse nonempty array of double real or double complex');
      end
      if (ijb==2)
	if (n~=size(A,1) || n~=size(A,2))
	  error('MATLAB:agmg:NotSameDim', ...
                'when IJOB==2 or IJOB==102, the input matrix A must have same dimensions as on previous call with IJOB==1 or IJOB==101');
        end
	if (isreal(A) && notcpl==0)
	  error('MATLAB:agmg:NotSameType', ...
                'when IJOB==2 or IJOB==102, the input matrix A must have same type (real or complex) as on previous call with IJOB==1 or IJOB==101');
	end
      else
        n = size(A,1);
        if (size(A,2) ~= n)
           error('MATLAB:agmg:NonSquareMatrix', 'the input matrix A must be square');
        end
        if(isreal(A)), notcpl=1; else, notcpl=0; end
        preprocessed=0;
      end
   end
   if (ijb == 0 || ijb >= 2)
      if ~isequal(size(b),[n,1])
          error('MATLAB:agmg:RSHsizeMatchCoeffMatrix', ...
                'the input right hand side B must be a column vector with the same number of rows as the input matrix A');
      end
      if issparse(b)
          error('MATLAB:agmg:RSHsizeMustBeFull', ...
          ['the input right hand side B cannot be a sparse vector' ...
          ' Please convert it [b=full(b);] before calling AGMG']);
      end
   end
   if (notcpl)
      if (ijb == 0 || ijb >= 2)
        if(~isreal(b))
	     error('MATLAB:agmg:bcomplexwhenAreal', ...
                   'for real input matrix A, the right hand side vector B must be real');
        end
        if (isempty(x0))
           [x,iter,resvec]=dmtlagmg(A,b,verbose,icg,maxit,tol,ijob );
        else
           if(ijb == 3)
               disp('warning: the input argument X0 is ignored when IJOB==3')
           elseif(issparse(x0))
               error('MATLAB:agmg:X0sizeMustBeFull', ...
                     'the input initial approximation X0 cannot be a sparse vector');
           elseif (~isequal(size(x0),[n,1]))
               error('MATLAB:agmg:X0sizeMatchCoeffMatrix', ...
                     'the input initial approximation X0 must be a column vector with the same number of rows as the input matrix A');
           elseif(~isreal(x0))
	       disp('Warning: the input matrix A and right hand side vector B are real;')
               disp('         only the real part of the initial guess X0 will be taken into account')
           end
           [x,iter,resvec]=dmtlagmg(A,b,verbose,icg,maxit,tol,ijob ,x0);
        end
      else
         dmtlagmg(A,b,verbose,icg,maxit,tol,ijob );
      end
   else
      if (~preprocessed)
         disp('Warning: complex version is tentative only')
      end
      if (ijb == 0 || ijb >= 2)
        if (isempty(x0))
           [x,iter,resvec]=zmtlagmg(A,b,verbose,icg,maxit,tol,ijob );
        else
           if(ijb == 3)
               disp('Warning: the input argument X0 is ignored when IJOB==3')
           elseif(issparse(x0))
               error('MATLAB:agmg:X0sizeMustBeFull', ...
                     'the input initial approximation X0 cannot be a sparse vector');
           elseif (~isequal(size(x0),[n,1]))
               error('MATLAB:agmg:X0sizeMatchCoeffMatrix', ...
                     'the input initial approximation X0 must be a column vector with the same number of rows as the input matrix A');
           end
           [x,iter,resvec]=zmtlagmg(A,b,verbose,icg,maxit,tol,ijob ,x0);
        end
      else
           zmtlagmg(A,b,verbose,icg,maxit,tol,ijob );
      end
   end
   flag=0;
   if (ijb == 1 || ijb == -1)
      iter=[];
      if (ijb == 1), preprocessed=1; else, preprocessed=0; end
      relres=[];
      x=[];
      resvec=[];
   elseif (ijb == 3)
      relres=[];
      resvec=[];
   else
      if (iter < 0), flag=1; iter=-iter; end
      if (nargout > 2), relres=norm(A*x-b)/norm(b); else, relres=[]; end
   end
end
