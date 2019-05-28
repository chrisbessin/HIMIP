function y = normalize(x)

n = size(x);
minX = repmat(min(x), [n(1), 1]);
maxX = repmat(max(x), [n(1), 1]);
y = (x - minX)./(maxX - minX);