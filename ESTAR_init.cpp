#include "ESTAR.h"

using namespace std;

ESTAR::ESTAR(int dimx,int dimy)
{
	Estar_init(&estar,dimx,dimy);
}

ESTAR::~ESTAR()
{
	Estar_fini(&estar);
}

void  ESTAR::Estar_init(estar_t *estar,size_t dimx,size_t dimy)
{
  
	std::fstream fs;
	fs.open("FILE/ESTAR_INIT.txt", std::fstream::in|std::fstream::out|std::fstream::trunc);
	if(fs.is_open())
	{
	  fs<< "INSIDE ESTAR_INIT function.\n" << endl;
	  fs.close();
	}
	else
	{
	  cout << "File not created";
	}
	//file1.open(filename,ios::in|ios::out);
	
	//file1 << "This is the Estar_init function." << endl;
	Estar_grid_init(&estar->grid,dimx,dimy);
	Estar_pqueue_init(&estar->pq,dimx+dimy);
	
	for(int i = 0; i < dimx; ++i )
	{
		for(int j = 0; j < dimy; ++j )
		{	
			Estar_set_speed(estar,i,j,1.0);
		}
	}
	Estar_set_goal(estar,GOALX,GOALY);
	have_goal = 1;
	//file1.close();
}



void ESTAR::Estar_pqueue_init(estar_pqueue_t * pq, size_t cap)
{
  	std::fstream file;
	file.open("FILE/Estar_pqueue_init.txt", std::fstream::in|std::fstream::out|std::fstream::trunc);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	
	pq->heap = (estar_cell_t**) malloc (sizeof(estar_cell_t**) * (cap+1));
	pq->len = 0;
	pq->cap = cap;
	file << "pq->heap = " << pq->heap << " pq->len = "<<pq->len << " pq->cap = " << pq->cap << endl;
	file.close();
}

void ESTAR::Estar_grid_init(estar_grid_t * grid, size_t dimx, size_t dimy)
{
  
	std::fstream file;
	file.open("FILE/Estar_grid_init.txt", std::fstream::in|std::fstream::out|std::fstream::trunc);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	
	size_t x,y;
	estar_cell_t *cell;
	estar_cell_t **nbor;
	size_t a = sizeof(estar_cell_t);
	grid->cell =(estar_cell_t *) malloc (a * dimx * dimy);
	
	grid->dimx = dimx;
        grid->dimy = dimy;
	file << "Grid->dimx = " << dimx << " Grid->dimy = " << dimy << " Grid->cell = " << grid->cell << endl; 
	int sum =0;
	for (x = 0; x < dimx; ++x)
        {
    		for (y = 0; y < dimy; ++y) 
		{
			//++sum;
			//cout << "Allocated memory " << sum << endl;
      			cell = estar_grid_at(grid, x, y);
			file << "X = " << x << " Y = " << y << " cell = " << cell << endl;
      			cell->cost = 1.0;
           		cell->phi = INFINITY;
            		cell->rhs = INFINITY;
            		cell->key = INFINITY;
            		cell->pqi = 0;
            		cell->flags = 0;

            		nbor = cell->nbor;
            		if (x > 0) 
	        	{/* west */
				*(nbor++) = cell - 1;
      			}
      			if (x < dimx - 1)
	  		{/* east */
				*(nbor++) = cell + 1;
      			}
      			if (y > 0)
	  		{/* south */
				*(nbor++) = cell - dimx;
      			}
      			if (y < dimy - 1)
	  		{/* north */
				*(nbor++) = cell + dimx;
      			}
      			*nbor = 0;
      
      			nbor = cell->prop;
      			if (x > 0) 
	  		{
				if (y > 0) 
				{/* south-west */
	  				*(nbor++) = cell - 1;
	  				*(nbor++) = cell - dimx;
				}
				if (y < dimy - 1)
				{/* north-west */
	  				*(nbor++) = cell - 1;
	  				*(nbor++) = cell + dimx;
				}
      			}
      			if (x < dimx - 1)
	  		{
				if (y > 0)
				{/* south-east */
					*(nbor++) = cell + 1;
	  				*(nbor++) = cell - dimx;
				}
				if (y < dimy - 1)
				{/* north-east */
	  				*(nbor++) = cell + 1;
	 				*(nbor++) = cell + dimx;
				}
      			}
      			*nbor = 0;
    		}
     	}
     	file.close();
}


void ESTAR::Estar_fini(estar_t * estar)
{
	 Estar_pqueue_fini (&estar->pq);
 	 Estar_grid_fini (&estar->grid);
}

void ESTAR::Estar_pqueue_fini (estar_pqueue_t * pq)
{
	cout << "Deallocating the memory of pqueue" << endl;
	free(pq->heap);
  	pq->len = 0;
  	pq->cap = 0;
}

void ESTAR::Estar_grid_fini (estar_grid_t * grid)
{
	cout << "Deallocating the memory of grids" << endl;
  	free(grid->cell);
  	grid->dimx = 0;
 	grid->dimy = 0;
}


void ESTAR::Estar_set_speed(estar_t *estar, size_t ix, size_t iy, double speed)
{
  	std::fstream file;
	file.open("FILE/Estar_set_speed.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	//cout << " Inside set speed \n"; 
	double cost;
  	estar_cell_t * cell;
  	estar_cell_t ** nbor;
  	cell = estar_grid_at (&estar->grid, ix, iy);
	file << " cell = " << cell ;
	// XXXX I'm undecided yet whether this check here makes the most
   	// sense. The other option is to make sure that the caller doesn't
        // place obstacles into a goal cell. The latter somehow makes more
 	// sense to me at the moment, so in gestar.c there is code to filter
  	// goal cells from the obstacle setting routines.
  	////  if (cell->flags & ESTAR_FLAG_GOAL) {
  	////    return;
  	////  }
	file << " cell->phi = " << cell->phi
	     << " cell->rhs = " << cell->rhs
	     << " cell->cost = " << cell->cost
	     << " cell->key = " << cell->key
	     << " cell->pqi = " << cell->pqi << endl;
  	if (speed <= 0.0) 
  	{
   		cost = INFINITY;
  	}
  	else 
	{
    		cost = 1.0 / speed;
  	}
  	if (cost == cell->cost) 
	{
    		return;
  	}
  
  	cell->cost = cost;
  	if (speed <= 0.0) 
	{
    		cell->phi = INFINITY;
    		cell->rhs = INFINITY;
    		cell->flags |= ESTAR_FLAG_OBSTACLE;
  	}
  	else 
	{
    		cell->flags &= ~ESTAR_FLAG_OBSTACLE;
  	}

	Estar_update(estar, cell);
  	for (nbor = cell->nbor; *nbor != 0; ++nbor) 
  	{
    		Estar_update (estar, *nbor);
  	}
  	file.close();
}

void ESTAR::Estar_set_goal (estar_t * estar, size_t ix, size_t iy)
{
    	std::fstream file;
	file.open("FILE/Estar_set_goal.txt", std::fstream::in|std::fstream::out|std::fstream::trunc);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	estar_cell_t * goal = estar_grid_at (&estar->grid, ix, iy);
	file << "Goal = " << goal << endl;
	goal->rhs = 0.0;
	goal->flags |= ESTAR_FLAG_GOAL;
	//printf("goal->flags = %d\n", goal->flags);
	goal->flags &= ~ESTAR_FLAG_OBSTACLE;
	file << "goal->flags = " << goal->flags << endl;
	//printf("goal->flags = %d\n", goal->flags);
	printf("Inserting the goal cell in the priority queue\n");
	Estar_pqueue_insert_or_update (&estar->pq, goal);
	
	file.close();
}


void ESTAR::Estar_update(estar_t * estar, estar_cell_t * cell)
{
	//cout << "Inside Estar_update.." << endl;
	// XXXX check whether obstacles actually can end up being
     	//updated. Possibly due to effects of estar_set_speed? 
     	std::fstream file;
	file.open("FILE/Estar_update.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	//file << " &estar->pq = " << &estar->pq << " cell = " << cell << endl;
  	if (cell->flags & ESTAR_FLAG_OBSTACLE) 
	{
		file << " &estar->pq = " << &estar->pq << " cell = " << cell << endl;
    		Estar_pqueue_remove_or_ignore (&estar->pq, cell);
    		return;
  	}
  
  	// Make sure that goal cells remain at their rhs, which is supposed
     	//to be fixed and only serve as source for propagation, never as
     	//sink. 
  	if ( !(cell->flags & ESTAR_FLAG_GOAL)) 
	{
		
    		calc_rhs (cell, Estar_pqueue_topkey (&estar->pq));
  	}
  
  	if (cell->phi != cell->rhs) 
	{
    		Estar_pqueue_insert_or_update (&estar->pq, cell);
  	}
  	else
	{
    		Estar_pqueue_remove_or_ignore (&estar->pq, cell);
 	}
 	
 	file.close();
}



void ESTAR::Estar_pqueue_insert_or_update (estar_pqueue_t * pq, estar_cell_t * cell)
{
       	std::fstream file;
	file.open("FILE/Estar_pqueue_insert_or_update.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
  	size_t len;
  	estar_cell_t ** heap;
	file << "pq->heap = " << pq->heap << endl;
  	//printf("cell->pqi = %zu\n",cell->pqi);
  	if (0 != cell->pqi) 
  	{
	      file << " cell->pqi = " << cell->pqi << " CALC_KEY = " << CALC_KEY(cell) << endl;
	      cell->key = CALC_KEY(cell);
	      // could probably make it more efficient by only bubbling down when
	      // the bubble up did not change cell->pqi
	      bubble_up (pq->heap, cell->pqi);
	      bubble_down (pq->heap, pq->len, cell->pqi);
	      return;
  	}
  
 	// grow heap, realloc if necessary
  	len = pq->len + 1;
  	if (len <= pq->cap) 
  	{
		//printf("pq->heap = %d\n",pq->heap);
    		heap = pq->heap;
  	}
  	else 
  	{
	    size_t cap;
	    cap = 2 * pq->cap;
	    heap =(estar_cell_t**) realloc (pq->heap, sizeof(estar_cell_t**) * (cap+1));
	    if (NULL == heap) 
	    {
      		errx (EXIT_FAILURE, __FILE__": %s: realloc", __func__);
	    }
	    pq->heap = heap;
	    pq->cap = cap;
  	}
  	pq->len = len;
  
 	// append cell to heap and bubble up
  	//printf("cell->rhs = %f,cell->phi = %f\n",cell->rhs,cell->phi);
  	cell->key = CALC_KEY(cell);
  	heap[len] = cell;
  	cell->pqi = len;		// initialize pqi 
  	bubble_up (heap, len);
	//printf("pqueue updated.\n");
	file.close();
}


void ESTAR::Estar_pqueue_remove_or_ignore (estar_pqueue_t * pq, estar_cell_t * cell)
{
	if (0 == cell->pqi) 
  	{
    		// This could be done by the caller for efficiency, but it is much
    		// more convenient to do it here.
   		 return;
  	}
  
  	pq->heap[cell->pqi] = pq->heap[pq->len];
  	pq->heap[cell->pqi]->pqi = cell->pqi; // keep pqi consistent! 
  	--pq->len;
  	bubble_down (pq->heap, pq->len, cell->pqi);
  	cell->pqi = 0;		// mark cell as not on queue 
}

void ESTAR::swap (estar_cell_t ** aa, estar_cell_t ** bb)
{
	  size_t ti;
  	  estar_cell_t *tc;
  	  ti = (*aa)->pqi;
  	  (*aa)->pqi = (*bb)->pqi;
      	  (*bb)->pqi = ti;
      	  tc = (*aa);
          (*aa) = (*bb);
          (*bb) = tc;
}

void ESTAR::bubble_up (estar_cell_t ** heap, size_t index)
{
	//cout << "Inside the Estar bubble_up function." << endl;
      	std::fstream file;
	file.open("FILE/Estar_bubble_up.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	size_t parent;
  	parent = index / 2;
  	//printf("parent = %zu\n",parent);
  	//printf("heap[index]->key = %zu\n",heap[index]->key);
  	//printf("heap[parent]->key = %zu\n",heap[parent]->key);	
  	while ((parent > 0) && (heap[index]->key < heap[parent]->key)) 
  	{
		//printf("into while of buble up\n");
    	swap (&heap[index], &heap[parent]);
    	index = parent;
    	parent = index / 2;
  	}
  	
  	file.close();
}

void ESTAR::bubble_down (estar_cell_t ** heap, size_t len, size_t index)
{
	//cout << "Inside the bubble_down function." <<  endl;
    	std::fstream file;
	file.open("FILE/Estar_bubble_down.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	size_t child, target;
  
  	target = index;
  	while (1)
	{
	  child = 2 * index;
	  file << " child = " << child << ",target = " << target << ",index = " << index << endl;
	  if (child <= len && heap[child]->key < heap[target]->key)
	  {
      		target = child;
	  }
	  ++child;
	  if (child <= len && heap[child]->key < heap[target]->key)
	  {
      		target = child;
	  }
	  if (index == target)
	  {
     		 break;
	  }
	  swap (&heap[target], &heap[index]);
	  index = target;
	}
	file.close();
}

double ESTAR::Estar_pqueue_topkey (estar_pqueue_t * pq)
{
	std::fstream file;
	file.open("FILE/Estar_pqueue_topkey.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	if (pq->len > 0) 
	{
		file << " pq->heap[1] = " << pq->heap[1]->key << endl; 
    		return pq->heap[1]->key;
  	}
  	return INFINITY;
	file.close();
}

void ESTAR::calc_rhs (estar_cell_t * cell, double phimax)
{
	std::fstream file;
	file.open("FILE/Estar_calc_rhs.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	estar_cell_t ** prop;
  	estar_cell_t * primary;
  	estar_cell_t * secondary;
  	double rr;
	
	file << " cell->rhs = " << cell->rhs << ", cell->prop = " << cell->prop << endl;
	
  	cell->rhs = INFINITY;
  	prop = cell->prop;
	file << "\t INSIDE WHILE" << endl;
  	while (NULL != *prop) 
	{
	//	file << "\t\t INSIDE WHILE"
    		primary = *(prop++);
		file << "\t\t primary = " << primary << " primary->rhs = " << primary->rhs << ", (*prop)->rhs) = " << (*prop)->rhs << endl;
		if (primary->rhs <= (*prop)->rhs)  
		{
			file << "\t\t Inside if statement " << " SEcondary = " << secondary << endl;
      			secondary = *(prop++);
    		}
    		else 
		{
			file << "\t\t Inside if statement " << endl;
      			secondary = primary;
      			primary = *(prop++);
    		}
    
    		// do not propagate from obstacles, queued cells, cells above the
    		// wavefront, or cells at infinity
    		if (primary->flags & ESTAR_FLAG_OBSTACLE
			|| primary->pqi != 0
			|| primary->phi > phimax
			|| isinf(primary->phi)) 
		{
     			 continue;
    		}
    
    		// the same goes from the secondary, but if that fails at least we
    		// can fall back to the non-interpolated update equation.
    		if (secondary->flags & ESTAR_FLAG_OBSTACLE
			|| secondary->pqi != 0
			|| secondary->phi > phimax
			|| isinf(secondary->phi))
		{
     		 	rr = primary->rhs + cell->cost;
    		}
    		else 
		{	
			file << "\t\t Going INSISE Estar_interpolate function.\n" << "\t\t cell->cost = "<< cell->cost
			                                                          <<  " primary->phi = " << primary-> phi
			                                                          <<  " secondary->phi = " << secondary->phi << endl;
      			rr = interpolate (cell->cost, primary->phi, secondary->phi);
			file << "\t\t\n rr = " << rr << endl;
    		}
                
    		if (rr < cell->rhs) 
		{
      			cell->rhs = rr;
    		}
  	}
  
  	if (isinf (cell->rhs)) 
	{
	    	// None of the above worked, we're probably done... but I have
    		// lingering doubts about about the effects of in-place primary /
    		// secondary sorting above, it could be imagined to create
    		// situations where we overlook something. So, just to be on the
    		// safe side, let's retry all non-interpolated options.
    		for (prop = cell->nbor; *prop != 0; ++prop) 
		{
      			rr = (*prop)->phi;
      			if (rr < cell->rhs)
			{
				cell->rhs = rr;
      			}
    		}
    		cell->rhs += cell->cost;
  	}
  	
  	file.close();
}

double ESTAR::interpolate (double cost, double primary, double secondary)
{
  	std::fstream file;
	file.open("FILE/Estar_interpolate.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	double tmp;
	file << " cost = " << cost << " primary = " << primary << " secondary = " << secondary << endl;
  	if (cost <= secondary - primary) 
	{
	    file << " Inside if statement : primary + cost = " << primary + cost << endl;
	    return primary + cost;
  	}
  
 	 // pow(cost,2) could be cached inside estar_set_speed. And so could
  	// the other squared terms. That might speed things up, but it would
  	// certainly make hearier caching code.
  
  	tmp = primary + secondary;
	file << " Return value  = " << (tmp + sqrt(pow(tmp, 2.0)
		                       - 2.0 * (pow(primary, 2.0)
			               + pow(secondary, 2.0)
			               - pow(cost, 2.0)))) / 2.0 << endl;
  	return (tmp + sqrt(pow(tmp, 2.0)
		     - 2.0 * (pow(primary, 2.0)
			      + pow(secondary, 2.0)
			      - pow(cost, 2.0)))) / 2.0;
			      
	file.close();
}

void ESTAR::Update()
{
     	std::fstream file;
	file.open("FILE/Estar_Update.txt", std::fstream::in|std::fstream::out|std::fstream::trunc);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	//cout << "Inside Update function." << endl;
	int status;
  	if (drag > 0 || estar.pq.len == 0) 
  	{
    		return;
  	} 
    	while (estar.pq.len != 0) 
	{ 
		//printf("estar.pq.len = %d.",estar.pq.len);
		file << "&estar = " << &estar << endl;
		Estar_propagate (&estar);
		if (estar.grid.dimx == STARTX && estar.grid.dimy == STARTY)
		{
			printf("Path obtained\n");
			break;
		}
  	}
  	status = Estar_check (&estar, "*** ");
  	if (0 != status) 
 	{
    		printf ("ERROR %d (see above)\n", status);
  	}
  	file.close();
}

void ESTAR::Estar_propagate (estar_t * estar)
{
     	std::fstream file;
	file.open("FILE/Estar_propogate.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	//cout << "Inside the Estar_propogate function." << endl;
	estar_cell_t * cell;
  	estar_cell_t ** nbor;
  
  	cell = Estar_pqueue_extract (&estar->pq);
	file << " cell = " << cell << endl;
  	if (NULL == cell) 
	{
    		return;
  	}
  
  	// The chunk below could be placed into a function called expand,
  	// but it is not needed anywhere else.
  	
  	if (cell->phi > cell->rhs) 
	{	
		file <<"Inside if statement\n" <<" cell->phi = " << cell->phi << " cell->rhs = " << cell->rhs << endl;
    		cell->phi = cell->rhs;
		file <<" cell->phi = " << cell->phi << " cell->rhs = " << cell->rhs << endl;
    		for(nbor = cell->nbor; *nbor != 0; ++nbor)
		{
			file << "\tgoing to update *nbor = " << *nbor << endl;
       			Estar_update (estar, *nbor);
    		}
  	}
  	else 
	{
		file <<"Inside else statement\n" <<" cell->phi = " << cell->phi << " cell->rhs = " << cell->rhs << endl;
    		cell->phi = INFINITY;
    		for (nbor = cell->nbor; *nbor != 0; ++nbor) 
		{
			file << "\tgoing to update *nbor = " << *nbor << endl;
      			Estar_update (estar, *nbor);
    		}
    		Estar_update (estar, cell);
  	}
  	file.close();
}


int  ESTAR::Estar_check (estar_t * estar, char const * pfx)
{
	int status = 0;
  	size_t ii, jj, kk;
  
  	for (ii = 0; ii < estar->grid.dimx; ++ii) 
	{
	    for (jj = 0; jj < estar->grid.dimy; ++jj) 
	    {
 	    	estar_cell_t * cell;
      		cell = estar_grid_at (&estar->grid, ii, jj);
      
      		if (cell->rhs == cell->phi) 
		{
			// consistent
			if (0 != cell->pqi) 
			{
	  			printf ("%sconsistent cell should not be on queue\n", pfx);
	  			status |= 1;
			}
      		}
      		else 
		{
			// inconsistent
			if (0 == cell->pqi) 
			{
	  			printf ("%sinconsistent cell should be on queue\n", pfx);
	  			status |= 2;
			}
      		}
      
      		if (0 == cell->pqi) 
		{
			// not on queue
			for (kk = 1; kk <= estar->pq.len; ++kk) 
			{
	  			if (cell == estar->pq.heap[kk]) 
				{
	    				printf ("%scell with pqi == 0 should not be on queue\n", pfx);
	    				status |= 4;
	    				break;
	  			}
			}
      		}
      		else 
		{
			// on queue
			for (kk = 1; kk <= estar->pq.len; ++kk) 
			{
	  			if (cell == estar->pq.heap[kk]) 
				{
	   				 break;
	  			}
			}
			if (kk > estar->pq.len) 
			{
	 			printf ("%scell with pqi != 0 should be on queue\n", pfx);
	  			status |= 8;
			}
      		}
	    }
	}
  
  	for (ii = 1; ii <= estar->pq.len; ++ii) 
	{
    		if (estar->pq.heap[ii]->pqi != ii) 
		{
      			printf ("%sinconsistent pqi\n", pfx);
      			Estar_dump_queue (estar, pfx);
      			status |= 16;
      			break;
    		}
  	} 
  	return status;
}

void ESTAR::Estar_dump_queue (estar_t * estar, char const * pfx)
{
	size_t ii;
  	for (ii = 1; ii <= estar->pq.len; ++ii) 
  	{
    		printf ("%s[%zu %zu]  pqi:  %zu  key: %g  phi: %g  rhs: %g\n",
	    			pfx,
	    			(estar->pq.heap[ii] - estar->grid.cell) % estar->grid.dimx,
	    			(estar->pq.heap[ii] - estar->grid.cell) / estar->grid.dimx,
	    			estar->pq.heap[ii]->pqi, estar->pq.heap[ii]->key,
	    			estar->pq.heap[ii]->phi, estar->pq.heap[ii]->rhs);
  	}
}


estar_cell_t* ESTAR::Estar_pqueue_extract (estar_pqueue_t * pq)
{
	//cout << "Inside the Estar_pqueue_extract" << endl;
	estar_cell_t * cell;
      	std::fstream file;
	file.open("FILE/Estar_pqueue_extract.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	} 
	file << " pq->len = " << pq->len; 
  	if (0 == pq->len) 
	{
    		return NULL;
  	}
  
 	cell = pq->heap[1];
	file << " cell = " << cell << endl;
  	cell->pqi = 0;		/* mark cell as not on queue */
  
  	if (1 == pq->len) 
	{
    		pq->len = 0;
    		return cell;
  	}
  
  	pq->heap[1] = pq->heap[pq->len];
  	pq->heap[1]->pqi = 1;		/* keep pqi consistent */
  	--pq->len;
  	// here would be a good place to shrink the heap
  
  	bubble_down (pq->heap, pq->len, 1);
	file.close();
  	return cell;
}

int ESTAR::Expose()
{
        std::fstream file;
	file.open("FILE/Estar_Expoase.txt", std::fstream::in|std::fstream::out|std::fstream::app);
	if (!file.is_open())
	{
	  cout << "Some error occured\n" << endl;
	}
	size_t ii, jj;
	double topkey, maxknown, maxoverall;
        estar_cell_t *cell;
	topkey = Estar_pqueue_topkey (&estar.pq);
  	maxknown = 0.0;
 	maxoverall = 0.0;
	printf("topkey = %f\n",topkey);
 	for (ii = 0; ii < DIMX; ++ii)
 	{
  	  for (jj = 0; jj < DIMY; ++jj)
  	  {
   	 	cell = estar_grid_at (&estar.grid, ii, jj);
		file << "cell = " << cell << " cell->phi = " << cell->phi << "cell->rhs = " << cell->rhs << endl; 
		 //printf("cell->pqi = %zu,cell->rhs = %zu\n",cell->pqi,cell->rhs);
  	        if (cell->rhs == cell->phi && isfinite(cell->rhs)) 
  	        {
		    ///file << "cell->pqi = " << cell->pqi <<",cell->rhs = " << cell->rhs << endl;
	            if (0 == cell->pqi && cell->rhs <= topkey && maxknown < cell->rhs) 
                    {
	                maxknown = cell->rhs;
	            }
	            if (maxoverall < cell->rhs)
                    {
	               maxoverall = cell->rhs;
	            }
	            file << " Maxknown = " << maxknown << ", maxoverall = " << maxoverall << endl;
                } 
          }
       }
       if (maxknown == 0.0) 
       {    
  	   maxknown = 0.0001;
       }
       if (maxoverall == 0.0) 
       {
           maxoverall = 0.0001;
       }
	printf("Maxknown = %.4f maxoverall = %.4f.\n",maxknown,maxoverall);
	/*If available, trace the path from start to goal*/
       cell = estar_grid_at (&estar.grid, STARTX, STARTY);
       file << " cell = " << cell << endl;
       file << "cell->pqi = " << cell->pqi <<" ,cell->rhs = " << cell->rhs  << endl;
       //printf("cell->pqi = %zu,cell->rhs = %zu\n",cell->pqi,cell->rhs);
       if (0 == cell->pqi && cell->rhs <= maxknown) 
       {
	    printf("Finding the traced path\n");
            double px, py, dd, dmax, ds;
            px = STARTX;
            py = STARTY;
            dmax = 1.3 * cell->rhs;
            printf ("dmax = %f,startX = %f,startY = %f\n",dmax,px,py);
            ds = 0.1;
            for (dd = 0.0; dd <= dmax; dd += ds) 
            {
                double gx, gy, gg;
                int ix, iy;
                if (0 == Estar_cell_calc_gradient (cell, &gx, &gy)) 
	        {
		    break;
                }
                gg = sqrt(pow(gx, 2.0) + pow(gy, 2.0));
	        //printf("gg = %f\n",gg);
      
	        gx *= ds / gg;
                gy *= ds / gg;
      
	        px += gx;
                py += gy;
		
	        ix = (int) rint(px);
                iy = (int) rint(py);
	        //printf("%d %d\n",ix,iy);
                if (ix < 0 || ix >= DIMX || iy < 0 || iy >= DIMY) 
	        {
		     break;
                }
                cell = estar_grid_at (&estar.grid, ix, iy);
      
	        if (cell->flags & ESTAR_FLAG_GOAL) 
	        {
		   break;
                }	
	   }
	   printf("Path is obtained.\n");	
  }
 
	return 1;
}

int ESTAR::Estar_cell_calc_gradient (estar_cell_t * cell, double * gx, double * gy)
{
	estar_cell_t ** nn;
  	estar_cell_t * n1;
  	estar_cell_t * n2;
  	int direction;
  
  	n1 = NULL;
  	for (nn = cell->nbor; *nn != NULL; ++nn) 
  	{
    		if (isfinite ((*nn)->rhs)
			&& (*nn)->rhs < cell->rhs
			&& (n1 == NULL || (*nn)->rhs < n1->rhs)) 
		{
      			n1 = *nn;
    		}
  	}
  	if (NULL == n1) 
	{
    		return 0;
  	}
  
  	direction = n1 - cell;
  	// +1 means right
  	// -1 means left
  	// +dimx means up (grid is arranged like pixels on a screen)
  	// -dimx means down
  
  	n2 = NULL;
  	for (nn = cell->nbor; *nn != NULL; ++nn) 
	{
    		if (isfinite ((*nn)->rhs)
			&& (*nn) != n1
			&& direction != cell - *nn /* check it is not opposite n1 */
			&& (n2 == NULL || (*nn)->rhs < n2->rhs)) 
		{
		      n2 = *nn;
    		}
  	}

 	if (NULL == n2) 
	{
    		if (direction == -1) 
		{
      			*gx = n1->rhs - cell->rhs; /* some negative value */
      			*gy = 0.0;
    		}
    		else if (direction == 1) 
		{
      			*gx = cell->rhs - n1->rhs; /* some positive value */
      			*gy = 0.0;
    		}
    		else if (direction < 0) 
		{
      			*gx = 0.0;
      			*gy = n1->rhs - cell->rhs; /* some negative value */
    		}
    		else 
		{
      			*gx = 0.0;
     		 	*gy = cell->rhs - n1->rhs; /* some positive value */
    		}
    		return 1;
  	}
  
  	if (direction == -1) 
	{
    		*gx = n1->rhs - cell->rhs;
    		if (cell - n2 > 0) 
		{
      			*gy = n2->rhs - cell->rhs;
    		}
    		else 
		{
      			*gy = cell->rhs - n2->rhs;
    		}	
  	}
  	else if (direction == 1) 
	{
    		*gx = cell->rhs - n1->rhs;
    		if (cell - n2 > 0) 
		{
      			*gy = n2->rhs - cell->rhs;
    		}
    		else 
		{
      			*gy = cell->rhs - n2->rhs;
    		}
  	}
  	else if (direction < 0) 
	{
    		if (cell - n2 > 0) 
		{
      			*gx = n2->rhs - cell->rhs;
    		}
    		else 
		{
      			*gx = cell->rhs - n2->rhs;
    		}
    		*gy = n1->rhs - cell->rhs;
  	}
  	else 
	{
    		if (cell - n2 > 0) 
		{
      			*gx = n2->rhs - cell->rhs;
    		}
    		else 
		{
      			*gx = cell->rhs - n2->rhs;
    		}
    		*gy = cell->rhs - n1->rhs;
  	}
  
  	return 2;
}

void ESTAR::Place_obstacle(int x,int y)
{
	mx = x;
	my = y;
	cout << "have_goal = " << have_goal << endl;
	if (mx >= 0 && mx < DIMX && my >= 0 && my < DIMY) 
	{
		estar_cell_t * cell = estar_grid_at(&estar.grid, mx, my);
    	        if (cell->flags & ESTAR_FLAG_GOAL) 
	        {
    	  		return;
   		}
		if (0 == have_goal) 
		{
      		    if (cell->flags & ESTAR_FLAG_OBSTACLE) 
	            {
				return;
      		    }
      		    Estar_set_goal (&estar, mx, my);
      		    have_goal = 1;
      			//gtk_widget_queue_draw (w_phi);
      		    return;
    	        }
		if (cell->flags & ESTAR_FLAG_OBSTACLE) 
		{
      		    drag = -1;
      		    change_obstacle (mx, my, ODIST, 0);
    	        }
    	        else
		{
      		    drag = -2;
      		    change_obstacle (mx, my, ODIST, 1);
                }
	}	
	return;
}

void ESTAR::change_obstacle (int cx, int cy, int dist, int add)
{
	int const dim = 2 * dist + 1;
  	double md2[dim * dim];
  	double *ptr;
  	double d2;
  	int x0, y0, x1, y1, x2, y2, x3, y3, ix, iy, jx, jy, nobst;
  
  	x0 = cx - dist;
  	if (0 > x0) 
  	{
   		x0 = 0;
  	}
  	y0 = cy - dist;
  	if (0 > y0) 
  	{
    		y0 = 0;
  	}
  	x1 = cx + dist + 1;
  	if (x1 > DIMX) 
  	{
    		x1 = DIMX;
  	}
  	y1 = cy + dist + 1;
  	if (y1 > DIMY) 
  	{
    		y1 = DIMY;
  	}
  
  	x2 = cx - 2 * dist;
  	if (0 > x2) 
  	{
    		x2 = 0;
  	}
  	y2 = cy - 2 * dist;
  	if (0 > y2) 
  	{
    		y2 = 0;
  	}
  	x3 = cx + 2 * dist + 1;
  	if (x3 > DIMX) 
  	{
    		x3 = DIMX;
  	}
  	y3 = cy + 2 * dist + 1;
  	if (y3 > DIMY) 
 	{
    		y3 = DIMY;
 	}
  
  	nobst = 0;
  	for (ix = x2; ix < x3; ++ix) 
  	{
    		for (iy = y2; iy < y3; ++iy) 
		{
      			if (0 == add && ix == cx && iy == cy) 
	  		{
				continue;
      			}
      			if ((0 != add && ix == cx && iy == cy)||(estar_grid_at (&estar.grid, ix, iy)->flags & ESTAR_FLAG_OBSTACLE))
	 		{
	  			ptr = md2;
	  			for (jx = x0; jx < x1; ++jx) 
				{
	    				for (jy = y0; jy < y1; ++jy) 
					{
	      					d2 = pow (ix-jx, 2.0) + pow (iy-jy, 2.0);
	      					if (0 == nobst || d2 < *ptr) 
						{
							*ptr = d2;
	      					}
	      					++ptr;
	    				}
	  			}
	  			++nobst;
			}
   			else 
			{
      
			}
    		}
  	}
  
  	if (0 == nobst)
  	{
    		for (ix = x0; ix < x1; ++ix) 
		{
      			for (iy = y0; iy < y1; ++iy)
	  		{
	  	 		Estar_set_speed (&estar, ix, iy, 1.0);
      			}
		}
  	}
  	else 
  	{
     		ptr = md2;
     		for (ix = x0; ix < x1; ++ix)
	 	{
        		for (iy = y0; iy < y1; ++iy)
			{
	       			d2 = sqrt(*(ptr++)) - 0.5;
	       			if (d2 < 0) 
		   		{
  			  		Estar_set_speed (&estar, ix, iy, 0.0);
	       			}
	       			else if (d2 >= dist) 
		   		{
	          			Estar_set_speed (&estar, ix, iy, 1.0);
	       			}
	       			else
		   		{
	          			Estar_set_speed (&estar, ix, iy, d2 / dist);
	       			}
      			}
    		}		
  	}
}
