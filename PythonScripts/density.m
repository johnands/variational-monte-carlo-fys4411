A = load('positionsNoInteraction.dat');

positionsXY = A(:,1:2);

colormap('gray');
[z,N] = hist3(positionsXY,[800 800]);
% surf(x)

xb = linspace(min(positionsXY(:,1)), max(positionsXY(:,1)), size(z,1)+1);
yb = linspace(min(positionsXY(:,2)), max(positionsXY(:,2)), size(z,1)+1);

imagesc(xb,yb,max(max(z))-z);
xlabel('x');
ylabel('y');
