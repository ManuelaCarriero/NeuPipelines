hf = figure('Units','normalized'); 
colormap gray
hCB = colorbar('north');
caxis([-1, 0]);
set(gca,'Visible',false)
hCB.Position = [0.15 0.3 0.74 0.4];
hf.Position(4) = 0.1000;