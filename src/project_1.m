%N=number of grid sections, there are N+1 nodes
%h is the spacing between nodes
N=30; h=(2*pi)/N;
x=0:h:2*pi;
y=0:h:2*pi;
u=ones(N+1,N+1);
%Values for F(x,y)
F=zeros(N+1,N+1);
for i=1:N+1
    for j=1:N+1
        F(i,j)=cos(x(i)/2+pi/2)*sin(y(j)/2);
    end
end
%Dirichlet Boundary Conditions
for i=1:1:N+1
    u(1,i)=(y(i))^3;
    u(N+1,i)=((y(i))^2)*cos(y(i));
    u(i,N+1)=((2*pi)^3)+x(i)*(2*pi-(2*pi)^2);
end
%5-Point Stencil Method (Gauss-Seidel)
%establishing starting values for error and # of iterations
h2=h^2;
err = 1;
iter=1;
%Error desired to be very small because of inaccuracy of L-infinity Norm
%While error is above threshold, program continues to iterate
while err>0.0000000001
    umaxprev=max(max(u(2:N,2:N)));
    iter=iter+1;
    %x varies from 2 to N as values at 1 and N+1 are fixed
    for i=2:N
        for j=2:N
            u(i,j)=(-F(i,j)*h2-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)))/-4;
        end
        %Neumann B.C.
        %To acheive du/dy=0, ghost node at node -1 is the same value as
        %node 2
        j=1;
        u(i,j)=(-F(i,j)*h2-(u(i+1,j)+u(i-1,j)+2*u(i,j+1)))/-4;
    end
    %Error specified with L-infinity Norm
    err = (max(max(u(2:N,2:N))-umaxprev))/umaxprev;
end
%Plotting the Results
surf(x,y,u)
xlabel('X'),ylabel('Y'),zlabel('U(X,Y)'), title('Gauss-Seidel Solution to Poisson Equation')