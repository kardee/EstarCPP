#include <iostream>
#include "ESTAR.h"

using namespace std;

int main()
{
	ESTAR E(DIMX,DIMY);
 	mx = -1;
  	my = -1;
  	drag = 0;
	
	E.Update();
	E.Expose();
	have_goal = 1;
	E.Place_obstacle(31,17);
	//cout << "END : " << endl;
	//E.Place_obstacle(46,24);
	E.Update();
	E.Expose();
	//cout << "END : " << endl;
	//E.Estar_fini(&estar);
	return 1;
}
