A = load('positionsInteraction6.dat');
%A1 = load('positionsNoInteraction6.dat');

positionsYZ = A(:,2:3);

pathname = '/home/fenics/Documents/FYS4411/variational-monte-carlo-fys4411/Report/';

colormap('bone');
[x,N] = hist3(positionsYZ,[800 800]);

yb = linspace(min(positionsYZ(:,1)), max(positionsYZ(:,1)), size(x,1)+1);
zb = linspace(min(positionsYZ(:,2)), max(positionsYZ(:,2)), size(x,1)+1);

figure(1)
imagesc(yb,zb,max(max(x))-x);
xlabel('z');
ylabel('y');
hold('on');
x1 = get(gca, 'xlim');
y1 = get(gca, 'ylim');
plot([-1,-1], y1, 'r-')
plot([1,1,], y1, 'r-')
plot(x1, [-1,-1], 'r-')
plot(x1, [1,1], 'r-')   
axis('equal');
title('Distribution of particles in yz-plane')
remove_frame(); 
print(gcf, '-dpdf', strcat(pathname, 'radialDist2YZ'));
