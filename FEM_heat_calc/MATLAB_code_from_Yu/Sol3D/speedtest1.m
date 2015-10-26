X=50; Y=50; Z=25;
T=1000;

%%%%% Flag
F=rand(X,Y,Z);
F(1,:,:)=2;
F(X,:,:)=2;
F(:,1,:)=2;
F(:,Y,:)=2;
F(:,:,1)=2;
F(:,:,Z)=2;


%%%%% For Loop
A=rand(X,Y,Z);

%%%%% No Loop
B=A;

%%%%% For Loop
tic
for t=1:T
    for x=1:X
        for y=1:Y
            for z=1:Z
                if (0<=F(x,y,z))&&(F(x,y,z)<0.5)
                    A(x,y,z)=A(x-1,y,z)+A(x+1,y,z)+A(x,y-1,z)+A(x,y+1,z)+A(x,y,z-1)+A(x,y+1,z+1);
                end
                if (0.5<=F(x,y,z))&&(F(x,y,z)<1)
                    A(x,y,z)=(A(x-1,y,z)+A(x+1,y,z)+A(x,y-1,z)+A(x,y+1,z)+A(x,y,z-1)+A(x,y+1,z+1))/2;
                end
            end
        end
    end
end
toc