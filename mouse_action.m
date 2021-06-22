function main
clc,close all 
f=figure;
aH=axes('Xlim',[0,1],'Ylim',[0,1]);
h=line([0.5 0.5], [0 1],'ButtonDownFcn',@startDragFcn);
set(f,'WindowButtonUpFcn',@stopDragFcn);


%-----------  vidines funkcijos ------------------
function startDragFcn(varargin)
set(f, 'WindowButtonMotionFcn',@draggingFcn)
end

function draggingFcn(varargin)
pt=get(aH,'Currentpoint');
set(h,'xData',pt(1)*[1 1]);
end

function stopDragFcn(varargin)
set(f, 'WindowButtonMotionFcn','')
end


end

%%%%%%%%%%%%%%%%%