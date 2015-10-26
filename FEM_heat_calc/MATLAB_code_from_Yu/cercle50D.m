function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 1.5788242632153422 1]);
set(ax,'PlotBoxAspectRatio',[2.5 0.63338271604938334 5]);
set(ax,'XLimMode','auto');
set(ax,'YLimMode','auto');
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pderect([-0.5 0 0.20000000000000001 -0.20000000000000001],'R1');
pdeellip(-0.25,0,0.050000000000000003,0.050000000000000003,...
0,'E1');
pderect([0 0.5 0.20000000000000001 -0.20000000000000001],'R2');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1-E1+R2')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(11,...
'dir',...
1,...
'1',...
'50')
pdesetbd(10,...
'dir',...
1,...
'1',...
'50')
pdesetbd(9,...
'dir',...
1,...
'1',...
'50')
pdesetbd(8,...
'dir',...
1,...
'1',...
'50')
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
'neu',...
1,...
'0',...
'0')
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
pdetool('refine')

% PDE coefficients:
pdeseteq(2,...
'1.0',...
'0.0',...
'0.0',...
'1.0',...
'0:2',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0';...
'0.0';...
'0.0';...
'1.0'])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
str2mat('0','1308','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 3 1 1 1 1 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');
