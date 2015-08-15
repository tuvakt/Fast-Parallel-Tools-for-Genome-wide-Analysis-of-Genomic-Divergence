/*
 * comparative.h
 *
 * Author: tuvakt
 */

#ifndef COMPARATIVE_H_
#define COMPARATIVE_H_

#include <sys/time.h>

int get_population_size(int *pos);
void slide_right(int *idx, int *positions, int start, int stop, int length);
int readln (int fd, char *buf, int max);
double time_ddiff(struct timeval before, struct timeval after);


#endif /* COMPARATIVE_H_ */
