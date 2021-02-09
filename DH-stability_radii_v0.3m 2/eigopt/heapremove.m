function heap = heapremove(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)


index = find(heap == index);

vertices(heap(index)).quad = inf;

temp = heap(index);
heap(index) = heap(heaplength);
heap(heaplength) = temp;

heaplength = heaplength - 1;




while ((2*index <= heaplength)		&	(vertices(heap(index)).quad > min(vertices(heap(2*index+1)).quad, vertices(heap(2*index)).quad)))

	if (vertices(heap(2*index+1)).quad < vertices(heap(2*index)).quad)
		temp = heap(2*index + 1);
		heap(2*index + 1) = heap(index);
		heap(index) = temp;

		index = 2*index+1;
	else
		temp = heap(2*index);
		heap(2*index) = heap(index);
		heap(index) = temp;

		index = 2*index;
	end
		
end


return;