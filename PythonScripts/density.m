%A = load('positionsInteraction6.dat');
%A1 = load('positionsNoInteraction6.dat');

positionsXY = A(:,2:3);

pathname = '/home/fenics/Documents/FYS4411/variational-monte-carlo-fys4411/Report/';

colormap('bone');
[z,N] = hist3(positionsXY,[800 800]);

xb = linspace(min(positionsXY(:,1)), max(positionsXY(:,1)), size(z,1)+1);
yb = linspace(min(positionsXY(:,2)), max(positionsXY(:,2)), size(z,1)+1);

figure(1)
imagesc(xb,yb,max(max(z))-z);
xlabel('x');
ylabel('y');
title('Distribution of particles in xy-plane')
remove_frame(); 
%print(gcf, '-dpdf', strcat(pathname, 'radialDist2XY'));
