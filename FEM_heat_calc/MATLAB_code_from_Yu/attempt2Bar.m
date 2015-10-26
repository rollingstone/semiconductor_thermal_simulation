function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[1.5 1 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-1 1]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pderect([-0.5 0 0.20000000000000001 -0.20000000000000001],'R1');
pderect([0 0.5 0.20000000000000001 -0.20000000000000001],'R2');
pderect([-0.5 0 0.0050000000000000001 -0.0050000000000000001],'R3');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-R3+R2')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(12,...
'neu',...
1,...
'0',...
'0')
pdesetbd(11,...
'neu',...
1,...
'0',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'20+(1-exp(-0.43*pi^2*t))/(1-exp(-0.43*pi^2*0.05))*10*0.05')
pdesetbd(7,...
'neu',...
1,...
'0',...
'0')
pdesetbd(6,...
'neu',...
1,...
'0',...
'0')
pdesetbd(5,...
'neu',...
1,...
'0',...
'0')
pdesetbd(4,...
'neu',...
1,...
'0',...
'0')
pdesetbd(3,...
'dir',...
1,...
'1',...
'20+(1-exp(-0.43*pi^2*t))/(1-exp(-0.43*pi^2*0.05))*10*0.05')
pdesetbd(2,...
'dir',...
1,...
'1',...
'20+(1-exp(-0.43*pi^2*t))/(1-exp(-0.43*pi^2*0.05))*10*0.05')
pdesetbd(1,...
'dir',...
1,...
'1',...
'20')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')

% PDE coefficients:
pdeseteq(2,...
'(x<=0).*0.43+(x>0).*5',...
'0.0',...
'0',...
'1.0',...
'0:3',...
'20.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['(x<=0).*0.43+(x>0).*5';...
'0.0                  ';...
'0                    ';...
'1.0                  '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
str2mat('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 1 0 4 1 1 1 1 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')
