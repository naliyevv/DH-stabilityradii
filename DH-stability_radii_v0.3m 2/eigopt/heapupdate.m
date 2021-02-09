function [heap, heaplength] = heapupdate(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)


% first find the vertex on the heap
index = find(heap == index);




while ((2*index <= heaplength)		&	(vertices(heap(index)).quad > min(vertices(heap(min(2*index+1,heaplength))).quad, vertices(heap(2*index)).quad)))

	if (vertices(heap(min(2*index+1,heaplength))).quad < vertices(heap(2*index)).quad)
		temp = heap(min(2*index + 1,heaplength));
		heap(min(2*index + 1,heaplength)) = heap(index);
		heap(index) = temp;

		index = min(2*index+1,heaplength);
	else
		temp = heap(2*index);
		heap(2*index) = heap(index);
		heap(index) = temp;

		index = 2*index;
	end
		
end


return;