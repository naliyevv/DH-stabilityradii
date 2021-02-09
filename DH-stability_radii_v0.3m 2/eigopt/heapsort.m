function heap = heapsort(vertices)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified on July 24th, 2012)



heap(1) = 1;
l = length(vertices);

for j = 2:l
	heap = heapinsert(heap,j-1,vertices,j);
end


return;