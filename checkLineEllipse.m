p = -1e-5;
q = 1e-5;

h = 0;
k = 0;
a = 3;
b = 1;
alpha = 0;

[x1 x2 y1 y2] = lineEllipseMatrix(a,b,h,k, alpha, p, q);

disp([x1 x2 y1 y2]);
disp(sqrt((x1-x2)^2 + (y1-y2).^2));
