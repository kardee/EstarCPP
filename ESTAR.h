#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <err.h>
#include <fstream>

#define DIMX 50
#define DIMY 50
#define GOALX 46
#define GOALY 23
#define ODIST 3
#define STARTX 19
#define STARTY 11
#include <fstream>
static int drag,have_goal,mx,my;
#define estar_grid_at(grid,ix,iy) (&(grid)->cell[(ix)+(iy)*(grid)->dimx])
#define CALC_KEY(cell) ((cell)->rhs < (cell)->phi ? (cell)->rhs : (cell)->phi)
using namespace std;

enum 
{
  ESTAR_FLAG_GOAL     = 1,
  ESTAR_FLAG_OBSTACLE = 2
};

typedef struct estar_cell_s 
{
  double cost;			 /* set this to 1/speed for "sensible" values */
  double phi;
  double rhs;
  double key;			 /* managed by pqueue */
  size_t pqi;			 /* managed by pqueue; pqi==0 means "not on queue" */
  int flags;
  struct estar_cell_s * nbor[5]; /* null-terminated array of neighbors */
  struct estar_cell_s * prop[9]; /* null-terminated array of pairwise propagators */
} estar_cell_t;

typedef struct 
{
  estar_cell_t * cell;
  size_t dimx, dimy;
} estar_grid_t;

typedef struct 
{
  estar_cell_t **heap;
  size_t len, cap;
} estar_pqueue_t;

typedef struct 
{
  estar_grid_t grid;   // estar_grid_t is structure included in grid.h
  estar_pqueue_t pq;   // estar_pqueue_t is structure included in pqueue.h
} estar_t;

#ifndef ESTAR_H
#define ESTAR_H


class ESTAR
{
	public:
		ESTAR(){};
		ESTAR(int dimx,int dimy);
		void Estar_init(estar_t *estar,size_t dimx,size_t dimy);
		void Estar_fini (estar_t * estar);
		void Update();
		int  Expose();
		void  Place_obstacle(int x,int y);
		~ESTAR();
		
	private:
		void Estar_pqueue_init(estar_pqueue_t * pq, size_t cap);
		void Estar_grid_init(estar_grid_t * grid, size_t dimx, size_t dimy);	
		void Estar_pqueue_fini (estar_pqueue_t * pq);
		void Estar_grid_fini (estar_grid_t * grid);
		void Estar_set_speed (estar_t *estar, size_t ix, size_t iy, double speed);
		void Estar_update(estar_t * estar, estar_cell_t * cell);
		void Estar_pqueue_insert_or_update (estar_pqueue_t * pq, estar_cell_t * cell);
		void Estar_pqueue_remove_or_ignore (estar_pqueue_t * pq, estar_cell_t * cell);
		void Estar_propagate (estar_t * estar);
		void Estar_set_goal (estar_t * estar, size_t ix, size_t iy);
		void Estar_dump_queue (estar_t * estar, char const * pfx);
		int  Estar_check (estar_t * estar, char const * pfx);
		int  Estar_cell_calc_gradient (estar_cell_t * cell, double * gx, double * gy);
		void swap (estar_cell_t ** aa, estar_cell_t ** bb);
		void bubble_up (estar_cell_t ** heap, size_t index);
		void bubble_down (estar_cell_t ** heap, size_t len, size_t index);
		double Estar_pqueue_topkey (estar_pqueue_t * pq);
		void calc_rhs (estar_cell_t * cell, double phimax);
		double interpolate (double cost, double primary, double secondary);
		estar_cell_t * Estar_pqueue_extract (estar_pqueue_t * pq);
		void change_obstacle (int cx, int cy, int dist, int add);
		estar_t estar;
		fstream file1; 
		//static int drag;
};

#endif
