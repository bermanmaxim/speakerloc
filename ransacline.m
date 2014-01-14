function [ p, inliers ] = ransacline( X, N, t )
%RANSACLINE Find line fitting X using N RANSAC trials. t controls the
% distance for the points being considered as outliers.
% Written by Maxim Berman
best=0;
inliers = [];
for iter=1:N,
    ind = randsample(size(X,2), 2);
    in = lineptdist(X(:, ind), X, t);
    if numel(in) > best,
        inliers = in;
        best = numel(inliers);
    end
end
% fprintf('%i inliers found\n', numel(inliers));
p = polyfit(X(1, inliers), X(2, inliers), 1);
end

function inliers = lineptdist(L, X, t)
% Returns the inliers at given distance t from a line
    x1 = L(:,1);
    x2 = L(:,2);
    n = length(X);
    d = zeros(n, 1);
    for i = 1:n
        p3 = X(:,i);
        lambda = dot((x2 - x1), (x2-x3)) / dot( (x1-p2), (x1-p2) );
        d(i) = norm(lambda * p1 + (1-lambda) * p2 - p3);
    end
    inliers = find(abs(d) < t);
end