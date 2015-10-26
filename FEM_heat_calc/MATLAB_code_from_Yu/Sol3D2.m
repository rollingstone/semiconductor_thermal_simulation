D = 1;
L = [0,1,0,1,0,1];
T = 1;
Mx = 20;
My = 20;
Mz = 20;
N = 100;
P = 50;
ox = (L(2)-L(1))/Mx;         % le pas de x
x = L(1)+[0:Mx]*ox;
oy = (L(4)-L(3))/My;
y = L(3)+[0:My]*oy;
oz = (L(6)-L(5))/Mz;
z = L(5)+[0:Mz]*oz;
ot = T/N;
t = [0:N]*ot;
u_xyz0 = inline('0','x','y','z');



% Initialisation
for i=1:Mx+1
    j=1:My+1;
    h=1:Mz+1;
    u(i,j,h)=u_xyz0(x(i),y(j),z(h));
end

rx = D*ot/(ox*ox);
ry = D*ot/(oy*oy);
rz = D*ot/(oz*oz);
rx1 = 1+2*rx;
rx2 = 1-2*rx;
ry1 = 1+2*ry;
ry2 = 1-2*ry;
rz1 = 1+2*rz;
rz2 = 1-2*rz;

for j=1:Mx-1
    A(j,j)=ry1;
    if j>1
        A(j-1,j)=-ry;
        A(j,j-1)=-ry;
    end
end

for i=1:My-1
    B(i,i)=rx1;
    if i>1
        B(i-1,i)=-rx;
        B(i,i-1)=-rx;
    end
end

for h=1:Mz-1
    C(i,i)=rz1;
    if h>1
        C(h-1,h)=-rz;
        C(h,h-1)=-rz;
    end
end

 for k=1:N
     t=k*ot;
     u_1=u;
     u(:,:,1)=u(:,:,2);
     u(:,:,Mz+1)=u(:,:,Mz);
     u(:,1,:)=u(:,2,:);
     u(:,Mx+1,:)=u(:,Mx,:);
     u(1,:,:)=u(2,:,:);
     u(My+1,:,:)=u(My,:,:);

          for h=1:9
%                uyx=reshape(u(:,:,h),My+1,Mx+1);
%                uyx_1=uyx;
              if mod(k,2)==0
                for y=2:My
                    xx=2:Mx;
                    for i=2:My                                 % Du bout au top
                       jj=2:Mx;
                       bx(i,jj)=[ry*u(i,1,h),zeros(1,Mx-3),ry*u(i,Mx+1,h)]+rx*(u_1(i-1,jj,h)+u_1(i+1,jj,h))+rx2*u_1(i,jj,h);
                    end
                    u(y,xx,h)=linsolve(A, bx(y,xx)');
%                    u(y,xx,h)=u(y,xx);
                end
              else
                 for x=2:Mx
                     yy=2:My;            
                     for j=2:Mx
                         ii=2:My;                                % Du bout au top                 
                         by(ii,j)=[rx*u(1,j,h);zeros(My-3,1);rx*u(My+1,j,h)]+ry*(u_1(ii,j-1,h)+u_1(ii,j+1,h))+ry2*u_1(ii,j,h);
                     end
                     u(yy,x,h) = linsolve(B,by(yy,x));
%                     u(yy,x,h)=u(yy,x);
                 end
              end
          end


          for h=10:11
%               u=reshape(u(:,:,h),My+1,Mx+1);
%               uyx_1=uyx;
            if mod(k,2)==0
                for y=2:My
                    xx=2:Mx;
                    for i=2:9                                 
                       jj=2:Mx;
                       bx(i,jj)=[ry*u(i,1,h),zeros(1,Mx-3),ry*u(i,Mx+1,h)]+rx*(u_1(i-1,jj,h)+u_1(i+1,jj,h))+rx2*u_1(i,jj,h);
                    end
                    for i=10:11
                        jj=2:Mx;
                        bx(i,jj)=[ry*u(i,1,h),zeros(1,Mx-3),ry*u(i,Mx+1,h)]+rx*(u_1(i-1,jj,h)+u_1(i+1,jj,h))+rx2*u_1(i,jj,h)+P*ot;
                    end
                    for i=12:My
                        jj=2:Mx;
                        bx(i,jj)=[ry*u(i,1,h),zeros(1,Mx-3),ry*u(i,Mx+1,h)]+rx*(u_1(i-1,jj,h)+u_1(i+1,jj,h))+rx2*u_1(i,jj,h);
                    end
                    u(y,xx,h)=linsolve(A, bx(y,xx)');
 %                   u(y,xx,h)=u(y,xx);
                end

            else
                for x=2:Mx                                 
                    yy=2:My;
                    for j=2:Mx
                        ii=2:9;
                        by(ii,j)=[rx*u(1,j,h);zeros(9-2,1)]+ry*(u_1(ii,j-1,h)+u_1(ii,j+1,h))+ry2*u_1(ii,j,h);
                    end
                    for j=2:Mx
                        ii=10:11;
                        by(ii,j)=zeros(2,1)+ry*(u_1(ii,j-1,h)+u_1(ii,j+1,h))+ry2*u_1(ii,j,h)+P*ot;
                    end
                    for j=2:Mx
                        ii=12:My;
                        by(ii,j)=[zeros(My-12,1);rx*u(My+1,j,h)]+ry*(u_1(ii,j-1,h)+u_1(ii,j+1,h))+ry2*u_1(ii,j,h);
                    end
                        u(yy,x,h) = linsolve(B,by(yy,x));
%                        u(yy,x,h)=u(yy,x);
                end
            end
          end

          for h=12:Mz
%               u=reshape(u(:,:,h),My+1,Mx+1);
%               uyx_1=uyx;
              if mod(k,2)==0
                for y=2:My
                    xx=2:Mx;
                    for i=2:My                                 
                       jj=2:Mx;
                       bx(i,jj)=[ry*u(i,1,h),zeros(1,Mx-3),ry*u(i,Mx+1,h)]+rx*(u_1(i-1,jj,h)+u_1(i+1,jj,h))+rx2*u_1(i,jj,h);
                    end
                    u(y,xx,h)=linsolve(A, bx(y,xx)');
%                    u(y,xx,h)=u(y,xx);
                end
              else
                 for x=2:Mx
                     yy=2:My;            
                     for j=2:Mx
                         ii=2:My;                                                 
                         by(ii,j)=[rx*u(1,j,h);zeros(My-3,1);rx*u(My+1,j,h)]+ry*(u_1(ii,j-1,h)+u_1(ii,j+1,h))+ry2*u_1(ii,j,h);
                     end
                     u(yy,x,h) = linsolve(B,by(yy,x));
%                     u(yy,x,h)=u(yy,x);
                 end
              end
          end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %   for i=1:9
        %       if mod(k,2)==0
        %         for z=2:Mz
        %             xx=2:Mx;
        %             for h=2:Mz                                
        %                jj=2:Mx;
        %                bx(h,jj)=[rz*u(i,1,h),zeros(1,Mx-3),rz*u(i,Mx+1,h)]+rx*(u_1(i,jj,h-1)+u_1(i,jj,h+1))+rx2*u_1(i,jj,h);
        %             end
        %             u(i,xx,z)=linsolve(A, bx(z,xx)');
        %         end
        %       else
        %          for x=2:Mx
        %              zz=2:Mz;            
        %              for j=2:Mx
        %                  hh=2:Mz;                                                
        %                  bz(hh,j)=[rx*u(i,j,1);zeros(Mz-3,1);rx*u(i,j,Mz+1)]+rz*(u_1(i,j-1,hh)+u_1(i,j+1,hh))+rz2*u_1(i,j,hh);
        %              end
        %              u(i,x,zz) = linsolve(C,bz(zz,x));
        %          end
        %       end
        %   end
        %   
        %   for i=10:11
        %       if mod(k,2)==0
        %         for z=2:Mz
        %             xx=2:Mx;
        %             for h=2:9                                 
        %                jj=2:Mx;
        %                bx(h,jj)=[rz*u(i,1,h),zeros(1,Mx-3),rz*u(i,Mx+1,h)]+rx*(u_1(i,jj,h-1)+u_1(i,jj,h+1))+rx2*u_1(i,jj,h);
        %             end
        %             for h=10:11
        %                 jj=2:Mx;
        %                 bx(h,jj)=[rz*u(i,1,h),zeros(1,Mx-3),rz*u(i,Mx+1,h)]+rx*(u_1(i,jj,h-1)+u_1(i,jj,h+1))+rx2*u_1(i,jj,h)+P*ot;
        %             end
        %             for h=11:Mz
        %              	jj=2:Mx;
        %                 bx(h,jj)=[rz*u(i,1,h),zeros(1,Mx-3),rz*u(i,Mx+1,h)]+rx*(u_1(i,jj,h-1)+u_1(i,jj,h+1))+rx2*u_1(i,jj,h);
        %             end
        %             u(i,xx,z)=linsolve(A, bx(z,xx)');
        %         end
        %          
        %     else
        %         for x=2:Mx                                 
        %             zz=2:Mz;
        %             for j=2:Mx
        %                 hh=2:9;
        %                 by(hh,j)=[rx*u(i,j,1);zeros(9-2,1)]+rz*(u_1(i,j-1,hh)+u_1(i,j+1,hh))+ry2*u_1(i,j,hh);
        %             end
        %             for j=2:Mx
        %                 ii=10:11;
        %                 by(hh,j)=zeros(2,1)+rz*(u_1(i,j-1,hh)+u_1(i,j+1,hh))+ry2*u_1(i,j,hh)+P*ot;
        %             end
        %             for j=2:Mx
        %                 hh=11:Mz;
        %                 bz(hh,j)=[zeros(Mz-11,1);rx*u(i,j,Mz+1)]+ry*(u_1(i,j-1,hh)+u_1(i,j+1,hh))+ry2*u_1(i,j,hh);
        %             end
        %                 u(i,x,zz) = linsolve(C,bz(zz,x));
        %         end
        %     end
        %   end
        %   
        %   for i=12:My
        %       if mod(k,2)==0
        %         for z=2:Mz
        %             xx=2:Mx;
        %             for h=2:Mz                                
        %                jj=2:Mx;
        %                bx(h,jj)=[rz*u(i,1,h),zeros(1,Mx-3),rz*u(i,Mx+1,h)]+rx*(u_1(i,jj,h-1)+u_1(i,jj,h+1))+rx2*u_1(i,jj,h);
        %             end
        %             u(i,xx,z)=linsolve(A, bx(z,xx)');
        %         end
        %       else
        %          for x=2:Mx
        %              zz=2:Mz;            
        %              for j=2:Mx
        %                  hh=2:Mz;                                                
        %                  bz(hh,j)=[rx*u(i,j,1);zeros(Mz-3,1);rx*u(i,j,Mz+1)]+rz*(u_1(i,j-1,hh)+u_1(i,j+1,hh))+rz2*u_1(i,j,hh);
        %              end
        %              u(i,x,zz) = linsolve(C,bz(zz,x));
        %          end
        %       end
        %   end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%            for j=1:Mx
%                uyz=reshape(u(:,j,:),My+1,Mz+1);
%                uyz_1=uyz;
%                if mod(k,2)==0
%                 for y=2:My
%                     zz=2:Mz;
%                     for i=2:9                                 % Duyz bouyzt auyz top
%                        hh=2:Mz;
%                        bz(i,hh)=[ry*uyz(i,1),zeros(1,Mz-3),ry*uyz(i,Mz+1)]+rz*(uyz_1(i-1,hh)+uyz_1(i+1,hh))+rz2*uyz_1(i,hh);
%                     end
%                     for i=10:11
%                         hh=2:9;
%                         bz(i,hh)=[ry*uyz(i,1),zeros(1,9-2)]+rz*(uyz_1(i-1,hh)+uyz_1(i+1,hh))+rz2*uyz_1(i,hh);
%                         hh=10:11;
%                         bz(i,hh)=zeros(1,2)+rz*(uyz_1(i-1,hh)+uyz_1(i+1,hh))+rz2*uyz_1(i,hh)+P*ot;
%                         hh=12:Mz;
%                         bz(i,hh)=[zeros(1,Mz-12),ry*uyz(i,Mz+1)]+rz*(uyz_1(i-1,hh)+uyz_1(i+1,hh))+rz2*uyz_1(i,hh);
%                     end
%                     for i=12:My
%                         hh=2:Mz;
%                         bz(i,hh)=[ry*uyz(i,1),zeros(1,Mz-3),ry*uyz(i,Mz+1)]+rz*(uyz_1(i-1,hh)+uyz_1(i+1,hh))+rz2*uyz_1(i,hh);
%                     end
%                     uyz(y,zz)=linsolve(A, bz(y,zz)');
%                     u(y,j,zz)=uyz(y,zz);
%                 end
% 
%             else
%                 for z=2:Mz                                 % De gauyzche a droite
%                     yy=2:My;
%                     for h=2:Mz
%                         ii=2:9;
%                         by(ii,h)=[rz*uyz(1,h);zeros(9-2,1)]+ry*(uyz_1(ii,h-1)+uyz_1(ii,h+1))+ry2*uyz_1(ii,h);
%                     end
%                     for h=2:9
%                         ii=10:11;
%                         by(ii,h)=zeros(2,1)+ry*(uyz_1(ii,h-1)+uyz_1(ii,h+1))+ry2*uyz_1(ii,h);
%                     end
%                     for h=10:11
%                         ii=10:11;
%                         by(ii,h)=zeros(2,1)+ry*(uyz_1(ii,h-1)+uyz_1(ii,h+1))+ry2*uyz_1(ii,h)+P*ot;
%                     end
%                     for h=12:Mz
%                         ii=10:11;
%                         by(ii,h)=zeros(2,1)+ry*(uyz_1(ii,h-1)+uyz_1(ii,h+1))+ry2*uyz_1(ii,h);
%                     end
%                     for h=2:Mz
%                         ii=12:My;
%                         by(ii,h)=[zeros(My-12,1);rz*uyz(My+1,h)]+ry*(uyz_1(ii,h-1)+uyz_1(ii,h+1))+ry2*uyz_1(ii,h);
%                     end
%                         uyz(yy,z) = linsolve(B,by(yy,z));
%                         u(yy,j,z)=uyz(yy,z);
%                 end
%              end  
%           end
end
[x,y,z]=meshgrid(1:21,1:21,1:21);
    uyz=reshape(u(:,3,:),My+1,Mz+1);
    uyx=reshape(u(:,:,10),My+1,Mx+1);
    uxz=reshape(u(3,:,:),Mx+1,Mz+1);
    surf(uxz);
    xlabel('x');
    ylabel('y');
    zlabel('u');
