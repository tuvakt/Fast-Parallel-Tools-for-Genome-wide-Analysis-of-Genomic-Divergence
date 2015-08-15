/*
 * comparative.c
 *
 * Author: tuvakt
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include "comparative.h"

#define FALSE 0
#define TRUE 1

/**
 * Get the size of the population
 * @param pos: an array with the positions of the genome
 * @return the size of the population
 */
int get_population_size(int *pos) {
    int firstpos = pos[0];
    int nextpos = firstpos;
    int size = 0;
    while(nextpos == firstpos) {
        size++;
	nextpos = pos[size];
    }
    return size;
}

/**
 * Get the array indices of the next window.
 * Slide the window wstep to the right.
 *
 * @param idx array of size 2 with the previous left/right indices.
 * 		Holds new left/right indices after the function has returned.
 * 		The window lies in indices [left, right-1] in the array.
 *
 * @param positions the array with positions in the genome
 * @param start the start position of the window
 * @param end the end position of the window
 * @param length the length of the positions array
 */
void slide_right(int *idx, int *positions, int start, int stop, int length) {

    /* idx: holds previous left/right idx
       when returned, holds new left/right idx */
    int left, right;
    
    left = idx[0];
    right = idx[1];
    
    // find the new array index position of left and right idx-pointer
    while (left < length && positions[left] < start) {
        left++;
    }
    
    while(right < length && positions[right] <= stop) {
        right++;
    }
    
    // we need the startidx for a and b, and the number of positions in the window
    idx[0] = left;
    idx[1] = right;
    
}

/*
 * author: Einar Sorensen
 * Read a line from a file
 */
int readln (int fd, char *buf, int max)
{
    int ant=0;
    int r;
    
    while (ant < max) {
        if ((r = read(fd, &buf[ant], 1))!=1) {
	    if (r < 0)
	        printf(" Error reading character: errno: %d read %d chars \n", errno, r);
	    else
	        printf("END OF FILE \n");
	    return -1;
	}
	if (buf[ant]=='\n' || buf[ant]=='\r') {
	    if (ant > 0) {
	        buf[ant] = 0;
	        break;
	    }
	}
	else
    	    ant++;
    }
    return ant;
}

/**
 * Get the time elapsed between two struct timevals
 * @param before the first struct timeval
 * @param after the last struct timeval
 */
double time_ddiff(struct timeval before, struct timeval after) {
    double sec_diff, msec_diff;
    
    sec_diff = (double)(after.tv_sec - before.tv_sec);
    msec_diff = (double)(after.tv_usec - before.tv_usec);
    
    return sec_diff + msec_diff/1000000.0;
}


