function [A,w,U,V,objHistory] = OSBC(X, nCluster, lambda1, lambda2)
% input
%     X:       nSmp * nDim, each column is a discrete attribute
%     Yl:      nSmp * 1, labeled data, unlabeled data is nan
%     dist:    seperable bregman divergence
%     lambda:  regularization term
%     nCluster:cluster number
% output
%     M:       nSmp * nSmp, concensus matrix
%     objHistory
%     Y_km: nSmp * nCluster via kmeans
%     Y_sr:   nSmp * nCluster via spectral rotation
%
if ~exist('s_type', 'var')
    s_type = 'inner';
end

if ~exist('c_type', 'var')
    c_type = 'kmeans';
end

if ~exist('maxIter', 'var')
    maxIter = 100;
end

[nSmp, nDim] = size(X);

%***************************************************
% Optimization problem
%  min (w, A, U, V) loss(A, U, V) + lambda ||A - sum_i wi Ai||^2 + lambda w^T S w
%  s.t. w \in delta
%***************************************************

%***************************************************
% Compute S
%***************************************************


S = zeros(nDim, nDim);
switch s_type
	case 'inner'
		for i1 = 1:nDim
			Ai = bsxfun(@ne, X(:,i1), X(:,i1)');
			for i2 = 1:nDim
				Aj = bsxfun(@ne, X(:,i2), X(:,i2)');
				S(i1,i2) = sum(sum(Ai.*Aj))/nSmp^2;
			end
        end
    case 'corr'
	otherwise
		error('not supported yet');
end

%***************************************************
% Init w
%***************************************************
w = ones(nDim,1)/nDim;

%***************************************************
% Init A
%***************************************************
A = zeros(nSmp, nSmp);
for i = 1:nDim
	Ai = bsxfun(@eq, X(:, i), X(:, i)');
	A = w(i) * Ai;
end

for iter = 1:maxIter
	%***************************************************
	% Step 1, given A, w, get U, V
	%***************************************************
	switch c_type
		case 'kmeans'
			[y, V] = litekmeans(A, nCluster, 'Replicates', 10);
			Y = LabelFormat(y);
		case 'sc'

		case 'smf'
		otherwise
			error('not supported yet');
	end
	%***************************************************
	% Step 2, given A, U, V, get w
	% min (w)  lambda1 * ||A - wi Ai||^2 + lambda 2 w^ S w
	% s.t.    w in delta 
	%***************************************************
	p = zeros(nDim,1);
	q = zeros(nDim,1);
	for i=1:nDim
		Ai = bsxfun(@eq, X(:,i), X(:,i)');
		p(i) = sum(sum(Ai.^2)); 
		q(i) = sum(sum(A .* Ai));
	end
	S2 = lambda2/lambda1 * (S + diag(p));
	w = quadprog(S2, -q,[],[],ones(1,nDim), 1);


	%***************************************************
	% Step 3, given w U, V, get A
	%***************************************************
	switch c_type
		case 'kmeans'
			A_merge = zeros(nSmp, nSmp);
			for i = 1:nDim
				Ai = bsxfun(@eq, X(:, i), X(:, i)');
				A_merge = w(i) * Ai;
			end
			A = (Y * V + lambda1 * A_merge)/(1 + lambda1);
		case 'sc'

		case 'smf'
		otherwise
			error('not supported yet');
	end
end