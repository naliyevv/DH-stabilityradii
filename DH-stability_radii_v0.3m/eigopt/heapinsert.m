function heap = heapinsert(heap,heaplength,vertices,index)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)



l = heaplength;

l = l+1;
heap(l) = index;



while (floor(l/2) > 0) & vertices(heap(l)).quad < vertices(heap(floor(l/2))).quad

	temp = heap(floor(l/2));
	heap(floor(l/2)) = heap(l);
	heap(l) = temp;

	l = floor(l/2);
end


return;