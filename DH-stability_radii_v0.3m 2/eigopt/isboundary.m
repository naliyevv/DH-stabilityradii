function bool = isboundary(vertices,vertnum)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified version July 30, 2012)
%

nonboundary = find(vertices(vertnum).index > 0);


if length(nonboundary) > 1
	bool = 0;
else
	bool = 1;
end


return;