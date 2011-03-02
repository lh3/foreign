/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

double ran1()
{
	return drand48();
}               
void seedit(char *flag)
{
	srand48(time(0) ^ getpid());
}
int commandlineseed(char **seeds)
{
	unsigned short seedv[3], *seed48();
	seedv[0] = atoi( seeds[0] );
	seedv[1] = atoi( seeds[1] );
	seedv[2] = atoi( seeds[2] );
	printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2]);
	seed48(seedv);
	return(3);
}
