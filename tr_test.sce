A=read("C:/MinGW/msys/1.0/home/vtorri/gitroot/sw/f1.dat", -1, 2);
B=read("C:/MinGW/msys/1.0/home/vtorri/gitroot/sw/f2.dat", -1, 2);
//plot(A(:,1), B(:,2));
//plot(A(:,1), [A(:,2), B(:,2)]);
plot(A(:,1), A(:,2) - B(:,2));
