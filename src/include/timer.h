#ifndef CUSTOM_TIMER_H
#define CUSTOM_TIMER_H 1
#include <sys/time.h>

typedef struct timeval Timer;

void reset(struct timeval *timer) {
    gettimeofday(timer, 0);
}

long lap(struct timeval timer) {
    Timer now;
    gettimeofday(&now, 0);
    return (now.tv_sec - timer.tv_sec) * 1000000 + (now.tv_usec - timer.tv_usec);
}
#endif
