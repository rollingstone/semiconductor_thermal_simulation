X=50;
Y=50;
Z=25;
M=zeros(X,Y,Z);
S='str';
tic
for i=1:1000
for x=1:X
    for y=1:Y
        for z=1:Z
            if strcmp(S,'str')
                M(x,y,z)=X+Y+Z;
            else
                M(x,y,z)=X+Y-Z;
            end
        end
    end
end
end
toc


X=50;
Y=50;
Z=25;
M=zeros(X,Y,Z);
tic
for i=1:1000
for x=1:X
    for y=1:Y
        for z=1:Z
            if x<=y
                M(x,y,z)=X+Y+Z;
            elseif x==y+1
                M(x,y,z)=X+Y-Z;
                elseif x==y+2
                M(x,y,z)=X+Y-Z;
                elseif x==y+3
                M(x,y,z)=X+Y-Z;
                elseif x==y+4
                M(x,y,z)=X+Y-Z;
                elseif x==y+5
                M(x,y,z)=X+Y-Z;
                elseif x==y+6
                M(x,y,z)=X+Y-Z;
                elseif x==y+7
                M(x,y,z)=X+Y-Z;
                elseif x==y+8
                M(x,y,z)=X+Y-Z;
                elseif x==y+9
                M(x,y,z)=X+Y-Z;
                elseif x==y+10
                M(x,y,z)=X+Y-Z;
            end
        end
    end
end
end
toc