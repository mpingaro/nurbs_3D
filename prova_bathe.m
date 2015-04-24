syms t

X = [1 2*t; 0 1];
dx = det(X);

S = [5770*t^2 3846*t;3846*t 13.462*t];

t = X*S*X'/dx;