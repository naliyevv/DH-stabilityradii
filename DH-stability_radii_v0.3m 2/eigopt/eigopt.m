function [f, z, parsout] = eigopt(funname,bounds,pars)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim 
% (Modified version July 23, 2013)
%
% call: [f, z, parsout] = eigopt(funname,bounds,parsin)
%
% task: globally optimize an eigenvalue function lambda(omega) over a given box
%
% Input:
%		funname (string) - name of the function evaluating 
%						   lambda(omega) and gradient of lambda(omega) at a given omega
%		bounds (struct)	 - bounds of the box on which the eigenvalue to be optimized
%						   bounds.lb (dx1 real)	 --- lower bounds
%						   bounds.ub (dx1 real)  --- upper bounds
%		pars (struct)	 - parameters
%						   pars.gamma (real)	 --- a lower bound for the minimum eigenvalue
%													 of the Hessian of lambda(omega) for all
%													 omega in the box. (if to be maximized, same but 
%													 for -lambda(omega))
%						   pars.tol (real)		 --- the tolerance for the accuracy of the computed
%													 globally optimal value of lambda over the box.
%						   pars.itertol(integer) --- maximum number of quadratic functions allowed
%						   pars.minmax (integer) --- minimize lambda(omega) if minmax = 0
%													 otherwise maximize lambda(omega)
%						   pars.isplot(integer)	 --- plots the underlying graph if isplot ~= 0 and 
%													 the problem is 2-dimensional; otherwise do not plot
%						   pars.isprint(integer) --- prints the details if isprint ~= 0;
%													 otherwise do not print
%						   pars.iskeyboard(integer)--- stops and waits for user at every 11 iterations
%													   if iskeyboard ~= 0; otherwise do not interact
%						   Also should have fields for the other parameters to be passed to funname 
% Output:
%		f (real)	   - computed globally optimal value of lambda(omega) on the box; differs from the exact 
%						 solution by no more than tol
%		z (dx1 real)   - computed global optimizer
%		parsout(struct)- parsout.nfevals (integer)					--- total number of function evaluations
%						 parsout.lbound (real)						--- a lower (upper) bound for the globally
%																		minimal (maximal) value of lambda(omega)
%																		over the box such that |f - lbound| < tol
%						 parsout.nvertices (array of integers)		--- number of vertices at every iteration
%						 parsout.newvertexlist (array of integers) --- number of newly added vertices at 
%																	every iteration
%						 parsout.deadvertexlist (array of integers)--- number of dead vertices at every iteration
%						 parsout.cpulist (array of reals) --- cputime required at every iteration


% KEY DATA STRUCTURES
% vertices(j).
%	adjacency - array of indices to adjacent vertices
%	adjnum - no of adjacent integers
%	index - array of indices indicating active constraints
%	coor - vector of dim keeping the coordinates
%	quad - value of the largest quadratic function
%
% heap
%	keeps the indices of vertices sorted from the largest to the
%	smallest according to the value of the maximal quadratic function




lb = bounds.lb;
ub = bounds.ub;




	if (isfield(pars,'gamma'))
		gamma = pars.gamma;
	else
		gamma = 2;
	end

	if (isfield(pars,'tol'))
		tol = pars.tol; 
	else
		tol = 10^-4;
	end

	if (isfield(pars, 'itertol'))
		itertol = pars.itertol;
	else
		itertol = 2000;
	end
		
	if (isfield(pars, 'minmax'))
		minmax = pars.minmax;
	else
		minmax = 0;
	end
		
	if (isfield(pars, 'isplot'))
		isplot = pars.isplot;
	else
		isplot = 0;
	end

	if (isfield(pars, 'isprint'))
		isprint = pars.isprint;
	else
		isprint = 0;
	end

	if (isfield(pars,'iskeyboard'))
		iskeyboard = pars.iskeyboard;
	else
		iskeyboard = 0;
	end




dim = length(lb);
dim2 = length(ub);


if (dim ~= dim2)
	error('lengths of lb and ub must be the same');
end





t1 = cputime;

for j=1:dim
	quad(j,1) = (lb(j) + ub(j))/2;
end

[fmin(1),gmin(:,1),~,tsdp,tchol] = feval(funname,quad(:,1),pars);



if (minmax ~= 0)
	fmin(1) = -fmin(1);
	gmin(:,1) = -gmin(:,1);
end

ubound = fmin(1);
z = quad(:,1);







for j = 1:2^dim

	k = j-1;
	l = 1;

	while (k ~= 0)
		if (mod(k,2) == 0)
			vertices(j).coor(l,1) = lb(l);
			vertices(j).index(l) = -(2*l - 1);
		else
			vertices(j).coor(l,1) = ub(l);
			vertices(j).index(l) = -2*l;
		end
		k = floor(k/2);
		l = l+1;
	end


	for k = l:dim
		vertices(j).coor(k,1) = lb(k);
		vertices(j).index(k) = -(2*k - 1);
	end




	vertices(j).adjnum = 0;

	vertices(j).quad = evalq(vertices(j).coor, quad(:,1), fmin(1), gmin(:,1),gamma);
	vertices(j).index(dim+1) = 1;

end





for j = 1:2^dim
	for k = j+1:2^dim
		if (length(intersect(vertices(j).index,vertices(k).index)) == dim)
			adjnum = vertices(j).adjnum;
			vertices(j).adjacency(adjnum+1) = k;
			vertices(j).adjnum = vertices(j).adjnum + 1;

			adjnum = vertices(k).adjnum;
			vertices(k).adjacency(adjnum+1) = j;
			vertices(k).adjnum = vertices(k).adjnum + 1;
		end
	end
end





heap = heapsort(vertices);
heaplength = 2^dim;

lbound = vertices(heap(1)).quad;

t2 = cputime;








iternum = 1;
boundaryl = 0;


cpulist = [t2-t1];
deadvertexlist = [0];
newvertexlist = [heaplength];
nvertices = [heaplength];


while ((ubound - lbound > tol) & (iternum < itertol))



	if (isprint ~= 0)
		if (minmax == 0)
			fprintf('iter:%d lowbound:%.8f upbound:%.8f heaplength:%d\n',iternum,lbound,ubound,heaplength);
		else
			fprintf('iter:%d lowbound:%.8f upbound:%.8f\n',iternum,-ubound,-lbound);
		end
	end


	t1 = cputime;
	iternum = iternum + 1;


	% center for the new quadratic model
	quad(:,iternum)	=	vertices(heap(1)).coor;

	% function value and gradient at the center of the new model
	[fmin(iternum),gmin(:,iternum),~,tsdpt,tcholt] = feval(funname,quad(:,iternum),pars);

    tsdp = tsdp + tsdpt;
    tchol = tchol + tcholt;

	if (minmax ~= 0)
		fmin(iternum) = -fmin(iternum);
		gmin(iternum) = -gmin(iternum);
	end


	if (fmin(iternum) < ubound)
		ubound = fmin(iternum);
		z = quad(:,iternum);
	end






	% (1) DETERMINE DEAD VERTICES
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% some of these are interior dead
	% some are dead on the boundary
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	stack(1) = heap(1);
	stackl = 1;


	notboundaryl = 0;
	notboundarylist(1) = -1;
	boundaryl = 0;
	boundarylist(1) = -1;


	while (stackl > 0)


		vertex = stack(stackl);
		stackl = stackl - 1;

		adjnum = vertices(vertex).adjnum;
		boundary  = 0;



		% EXPAND TO THE ADJACENT VERTICES
		for j = 1:adjnum

			adjacent = vertices(vertex).adjacency(j);


			if (ismember(adjacent, union(union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl)),stack(1:stackl))) == 0)
			% otherwise adjacent is already inspected or to be inspected soon (in the stack)

				qnew = evalq(vertices(adjacent).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

				

				if (qnew > vertices(adjacent).quad)
					% ADJACENT IS DEAD, ADD ADJACENT TO THE STACK
					stackl = stackl + 1;
					stack(stackl) = adjacent;
				else
					% ADJACENT IS ALIVE, SO VERTEX IS ON THE BOUNDARY OF THE SET OF DEAD VERTICES
					boundary = 1;
				end
			end
		end





		if	boundary == 0
			% DEAD INTERIOR VERTICES
			notboundaryl = notboundaryl + 1;
			notboundarylist(notboundaryl) = vertex;
		else
			% DEAD BOUNDARY VERTICES
			boundaryl = boundaryl + 1;
			boundarylist(boundaryl) = vertex;
		end

	end









	deadboundaryl = 0;


	% (2) REMOVE INTERIOR DEAD VERTICES
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Easy, no need to modify the alive
	% vertices or introduce new vertices
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for j = 1:notboundaryl

		if isboundary(vertices,notboundarylist(j));
			deadisboundaryvertex = 1;
			deadboundaryl = deadboundaryl + 1;
			deadboundarylist(deadboundaryl) = notboundarylist(j);
		else
			deadisboundaryvertex = 0;
		end


		if (deadisboundaryvertex == 0)
			% Dead vertex is not a boundary vertex, remove it
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			heap = heapremove(heap,heaplength,vertices,notboundarylist(j));
			vertices(heap(heaplength)).quad = inf;

			heaplength = heaplength - 1;

		else
			% Dead vertex is a boundary vertex, cannot be removed
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% (i) but need to update the function value
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			deadboundary = notboundarylist(j);
			vertices(deadboundary).quad = evalq(vertices(deadboundary).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

			%%%%%%%%%%%%%%%%%%%%%%%%
			% (ii) and the index set
			%%%%%%%%%%%%%%%%%%%%%%%%
			for k = 1:dim+1
				if ( vertices(deadboundary).index(k) > 0 )
					vertices(deadboundary).index(k) = iternum;
					break;
				end
			end

			vertices(deadboundary).adjnum = 0;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% (iii) finally update the heap
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			heap = heapupdate(heap,heaplength,vertices,deadboundary);
		end
	end








	% (3) REMOVE DEAD VERTICES ON THE BOUNDARY
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Subtle, (i) introduce a new vertex in between
	% any pair of dead and alive vertices (ii) modify
	% the adjacency of an alive vertex next to a dead
	% vertex
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	verticesl = length(vertices);
	oldlength = verticesl;
	alldead = union(boundarylist(1:boundaryl),notboundarylist(1:notboundaryl));

	


	% GO THROUGH THE LIST OF DEAD VERTICES ON THE BOUNDARY
	for j = 1:boundaryl

		dead = boundarylist(j);
		adjnum = vertices(dead).adjnum;


		qnewdead = evalq(vertices(dead).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

		
		if isboundary(vertices,dead);
			deadisboundaryvertex = 1;
			deadboundaryl = deadboundaryl + 1;
			deadboundarylist(deadboundaryl) = dead;
		else
			deadisboundaryvertex = 0;
		end


		% Determine the adjacent alive vertices
		alivelist = setdiff(vertices(dead).adjacency, intersect(vertices(dead).adjacency, alldead));
		alivelength = length(alivelist);

		for l = 1:alivelength
			alive = alivelist(l);

			% THIS ADJACENT VERTEX IS ALIVE, MUST CREATE A NEW VERTEX


			qnewalive = evalq(vertices(alive).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum),gamma);

				

					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% I - FORM A NEW VERTEX BETWEEN THE DEAD and ALIVE
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					verticesl = verticesl + 1;
					%%%%%%%%%%%%%%%%%%%%
					%(i) its coordinates
					%%%%%%%%%%%%%%%%%%%%
					alpha = (qnewalive - vertices(alive).quad)/((qnewalive - vertices(alive).quad) - (qnewdead - vertices(dead).quad));
					vertices(verticesl).coor = alpha*vertices(dead).coor  +  (1-alpha)*vertices(alive).coor;


					%%%%%%%%%%%%%%%%%%%%
					%(ii) its index
					%%%%%%%%%%%%%%%%%%%%
					vertices(verticesl).index = intersect(vertices(alive).index,vertices(dead).index);
					vertices(verticesl).index(dim+1) = iternum;


					%%%%%%%%%%%%%%%%%%%%%%%%%
					%(iii) its function value
					%%%%%%%%%%%%%%%%%%%%%%%%%
					vertices(verticesl).quad = evalq(vertices(verticesl).coor, quad(:,iternum), fmin(iternum), gmin(:,iternum), gamma);
					for k = 1:iternum-1



						qval = evalq(vertices(verticesl).coor, quad(:,k), fmin(k), gmin(:,k), gamma);
						if (qval > vertices(verticesl).quad)
							vertices(verticesl).quad = qval;
						end
					end


					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					%(iv) initialize its adjacency list
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					vertices(verticesl).adjnum = 1;
					vertices(verticesl).adjacency(1) = alive;





					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% II - ADD THE NEW VERTEX TO THE HEAP
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					heap = heapinsert(heap,heaplength,vertices,verticesl);
					heaplength = heaplength + 1;

				



					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% III - UPDATE THE ADJACENCY LIST OF THE ALIVE VERTEX
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					k = 1;
					while (vertices(alive).adjacency(k) ~= dead)
						k = k+1;
					end
					vertices(alive).adjacency(k) = verticesl;


		end




		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% IV - REMOVE THE DEAD VERTEX
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if (deadisboundaryvertex == 0)
			% dead vertex is not a boundary vertex 
			heap = heapremove(heap,heaplength,vertices,dead);
			vertices(dead).quad = inf;

			heaplength = heaplength - 1;
		else
			% dead vertex is on the boundary
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% (i) need to update the function value
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			vertices(dead).quad = qnewdead;

			%%%%%%%%%%%%%%%%%%%%%%%%
			% (ii) and the index set
			%%%%%%%%%%%%%%%%%%%%%%%%
			for k = 1:dim+1
				if ( vertices(dead).index(k) > 0 )
					vertices(dead).index(k) = iternum;
					break;
				end
			end

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% (iii) and the adjacency
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			vertices(dead).adjnum = 0;

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% (iv) finally update the heap
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			heap = heapupdate(heap,heaplength,vertices,dead);
		end



	end







	% (4) FORM THE CONNECTIONS IN THE NEW POLYTOPE
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for j = oldlength+1:verticesl
		for k = j+1:verticesl

			if (vertices(j).adjnum == dim+1)
				break;
			end
		
			if (length(intersect(vertices(j).index,vertices(k).index)) == dim)
				% newly added jth and kth vertices are adjacent
				adjnum = vertices(j).adjnum;
				vertices(j).adjacency(adjnum+1) = k;
				vertices(j).adjnum = vertices(j).adjnum + 1;

				adjnum = vertices(k).adjnum;
				vertices(k).adjacency(adjnum+1) = j;
				vertices(k).adjnum = vertices(k).adjnum + 1;	
			end
		end
	end

	
	for j = oldlength+1:verticesl
		for k = 1:deadboundaryl
			if (length(intersect(vertices(j).index,vertices(deadboundarylist(k)).index)) == dim)
				% newly added jth and the kth boundary vertices are adjacent
				adjnum = vertices(j).adjnum;
				vertices(j).adjacency(adjnum+1) = deadboundarylist(k);
				vertices(j).adjnum = vertices(j).adjnum + 1;
		
				adjnum = vertices(deadboundarylist(k)).adjnum;
				vertices(deadboundarylist(k)).adjacency(adjnum+1) = j;
				vertices(deadboundarylist(k)).adjnum = vertices(deadboundarylist(k)).adjnum + 1;	
			end
		end
	end


	for j = 1:deadboundaryl
		for k = j+1:deadboundaryl
			if (length(intersect(vertices(deadboundarylist(j)).index,vertices(deadboundarylist(k)).index)) == dim)
			% dead jth and the kth boundary vertices are adjacent
				adjnum = vertices(deadboundarylist(j)).adjnum;
				vertices(deadboundarylist(j)).adjacency(adjnum+1) = deadboundarylist(k);
				vertices(deadboundarylist(j)).adjnum = vertices(deadboundarylist(j)).adjnum + 1;

				adjnum = vertices(deadboundarylist(k)).adjnum;
				vertices(deadboundarylist(k)).adjacency(adjnum+1) = deadboundarylist(j);
				vertices(deadboundarylist(k)).adjnum = vertices(deadboundarylist(k)).adjnum + 1;	
			end
		end
	end


	lbound = vertices(heap(1)).quad;





	t2 = cputime;

	cpulist = [cpulist; t2-t1];
	deadvertexlist = [deadvertexlist; boundaryl+notboundaryl];
	newvertexlist = [newvertexlist; verticesl-oldlength];
	nvertices = [nvertices; heaplength];



	% plot the graph (only in 2d) if asked for and every 10 iterations
	if (isplot & (dim == 2) & (mod(iternum,11) == 0))
		plot_dead(vertices,heap,heaplength,boundarylist,boundaryl,notboundarylist,notboundaryl,3);
	end

	if (isplot & (dim == 2) & (mod(iternum,10) == 0))
		display(iternum);
		plot_graph(vertices,heap,heaplength,boundarylist,boundaryl,notboundarylist,notboundaryl,3);
	end



	% interact with the user if asked for and every 10 iterations
	if (iskeyboard & (mod(iternum,11) == 0))
		keyboard;
	end
	

end

% keyboard


if (minmax == 0)
	f = ubound;
	parsout.lbound = lbound;
else
	f = -ubound;
	parsout.lbound = -lbound;
end


parsout.nfevals = iternum;
parsout.nvertices = nvertices;
parsout.deadvertexlist = deadvertexlist;
parsout.newvertexlist = newvertexlist;
parsout.cpulist = cpulist;
parsout.tsdp = tsdp;
parsout.tchol = tchol;


return;
