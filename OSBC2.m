% min ||A - \sum_k w_k A_k||^2 + 2*lambda*trace(F'* L_A *F)
%  st A>=0, A*1=1, F'*F=I, w>=0, 1^T w = 1  
% 
function [y, A, objHistory, lambdaHistory] = OSBC2(X, nCluster, varargin)

%Input 
%   X : nSmp * nDim, Rows of X correspond to points, 
%                    columns correspond to variables.
%   nCluster: number of cluster
%   [ ... ] = OSBC(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm
%   used by OSBC.  Parameters are:
%
%   'k' - number of neighbors to determine the initial graph, 
%         and the parameter r if r<=0
%         (default 15)
%   'r' - which could be set to a large enough value. 
%         If r<0, then it is determined by algorithm with k
%         (default -1)
%   'islocal' - determine the graph tobe local connected or all connected.  Choices are:
%            {1} - only update the similarities of the k neighbor pairs
%            {0} - update all the similarities
%         (default 1)
%   'maxiter' - Maximum number of iterations allowed.  Default is 100.
%
%   'y' - ground truth, only used for debug purpose
%         (default []) 

[nSmp, nDim] = size(X);

%**********************************************
% Parameters extraction
%**********************************************
pnames = {'k' 'r' 'islocal' 'maxiter'};
dflts =  {15  -1     1         50    };
if nargin < 2
    error('OSBC:TooFewInputs','At least two input arguments required.');
end
[eid,errmsg,initK,r,islocal,maxiter] = getargs(pnames, dflts, varargin{:});
if ~isempty(eid)
    error(sprintf('OSBC:%s',eid),errmsg);
end


%***************************************************
% Init w
%***************************************************
w = ones(nDim,1)/nDim;

%**********************************************
% Init A
%**********************************************
Aw = zeros(nSmp, nSmp);
for i = 1:nDim
	Ai = bsxfun(@eq, X(:, i), X(:, i)');
	Aw = w(i) * Ai + Aw;
end
Aw = Aw - diag(diag(Aw));
%**********************************************
% Optimization problem
% min ||a - Awi||^2
% s.t. a >=0 1' a = 1, nnz(a) = k
%**********************************************
if islocal
	for i = 1:nSmp
		A(i,:) = EProjSimplex_new(Aw(i,:), k);
	end
else
	for i = 1:nSmp
		A(i,:) = EProjSimplex_new(Aw(i,:));
	end	
end
%**********************************************
% Init Laplacian and F
%**********************************************
A0 = (A + A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[F, temp, evs] = eig1(L0, nCluster, 0);

if sum(evs(1:nCluster+1)) < 0.00000000001
    error('The original graph has more than %d connected component', nCluster);
end;

for iter = 1:maxiter
	%**********************************************
	% Update A (similarity matrix)
	%**********************************************
	Aw = zeros(nSmp, nSmp);
	for i = 1:nDim
		Ai = bsxfun(@eq, X(:, i), X(:, i)');
		Aw = w(i) * Ai + Aw;
	end
	ff = sum(F .* F, 2);
	distf = bsxfun(@plus, ff, ff') - 2 * (F * F');
	distf = max(real(distf), 0);
	%**********************************************
	% Optimization problem
	% min a^T a - 2 a^T Aw(i,:) + lambda a^T df
	% min ||a + s||^2
	% s = Awi - 0.5 * lambda * df
	%**********************************************
	S = Aw - 0.5 * lambda * distf;
	S = S - diag(diag(S));
	if islocal
		for i = 1:nSmp
			A(i,:) = EProjSimplex_new(S(i,:), k);
		end
	else
		for i = 1:nSmp
			A(i,:) = EProjSimplex_new(S(i,:));
		end	
	end
	%**********************************************
	% Update w (weight)
	%**********************************************
	p = zeros(nDim,1);
	q = zeros(nDim,1);
	for i=1:nDim
		Ai = bsxfun(@eq, X(:,i), X(:,i)');
		p(i) = sum(sum(Ai.^2));
		q(i) = sum(sum(A .* Ai));
	end
	w = quadprog(diag(p), -q, [], [], ones(1,m), 1);

	%**********************************************
	% Update Laplacian and F
	%**********************************************
    A = (A + A')/2;
    D = diag(sum(A));
    L = D - A;
    F_old = F;
    [F, temp, ev]=eig1(L, nCluster, 0);
    evs(:,iter+1) = ev;

	%**********************************************
	% Update regularization parameter
	%**********************************************
    fn1 = sum(ev(1:nCluster));
    fn2 = sum(ev(1:nCluster+1));
    if fn1 > 0.00000000001
    	% more than nCluster component
        lambda = 2*lambda;
    elseif fn2 < 0.00000000001
        lambda = lambda/2;  F = F_old;
        % less than nCluster compoent, back off
    else
        break;
    end;

end

%**********************************************
% Get discrete result
%**********************************************
[clusternum, y] = graphconncomp(sparse(A)); 
y = y';
if clusternum ~= nCluster
    sprintf('Can not find the correct cluster number: %d', nCluster)
end;

end

function [x ft] = EProjSimplex_new(v, k)
% 
%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=1
%  Reference: http://heimingx.cn/2016/10/15/euclidean-projection-onto-the-simplex/
%

if nargin < 2
    k = 1;
end;

ft=1;
n = length(v);

v0 = v-mean(v) + k/n;
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1>0;
        npos = sum(posidx);
        g = -npos;
        f = sum(v1(posidx)) - k;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0);

else
    x = v0;
end;
end