function quadval = plot_graph(vertices,heap,heaplength,boundarylist,boundaryl,notboundarylist,notboundaryl,fignum)



figure(fignum);
hold off;

for j = 1:heaplength
	
	cvertex = heap(j);
	
	for k = 1:vertices(cvertex).adjnum;
		avertex = vertices(cvertex).adjacency(k);	
		plot([vertices(cvertex).coor(1) vertices(avertex).coor(1)],[vertices(cvertex).coor(2) vertices(avertex).coor(2)],'k-');  

		if ((j == 1) & (k == 1))
			hold on;
		end
	end

	
	plot(vertices(cvertex).coor(1),vertices(cvertex).coor(2),'k.');
end



return;