function quadval = plot_dead(vertices,heap,heaplength,boundarylist,boundaryl,notboundarylist,notboundaryl,fignum)


figure(fignum);
hold on

for j = 1:boundaryl
	cvertex = boundarylist(j);
	plot(vertices(cvertex).coor(1),vertices(cvertex).coor(2),'rs');
end

for j = 1:notboundaryl
	cvertex = notboundarylist(j);
	plot(vertices(cvertex).coor(1),vertices(cvertex).coor(2),'rs');
end


 keyboard



return;