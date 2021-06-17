/*
-------------MIPT-------------
--Petrosyan Hakob, group 814--
--------No copyright.---------
*/

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>

/*
================================================================
= The program performs numerical integration on shared memory. =
=                   function: f(x) = cos(1/x)                  =
=                         0.01 < x < 10                        =
================================================================
*/

double LOWB = 0.01;
double HIGHB = 10.0;
int N_THREADS = 4;

double TRUE_VALUE = 8.47914;

double f(double x)
{
    return cos(1/x);
}

typedef struct {
    double invBegin_;
    double invEnd_;
    int nSteps_;
    double *bank_;
    pthread_mutex_t *mutex_;
    double timeElapsed_;
} thread_args_t;

typedef struct {
    thread_args_t args_;
    pthread_t pid_;
} thread_info_t;

void *compute(void *arg)
{
    struct timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    thread_args_t *threadArgs = arg;

    double lb = 1.0 / threadArgs->invEnd_;
    double hb = 1.0 / threadArgs->invBegin_;
    int n = threadArgs->nSteps_;
    double dx = (hb - lb) / n;

    double localResult = 0;
    for (int i = 0; i < n; ++i) {
        localResult += f(lb + i * dx) + f(lb + (i + 1) * dx);
    }
    localResult *= dx / 2;

    pthread_mutex_lock(threadArgs->mutex_);
    *(threadArgs->bank_) += localResult;
    pthread_mutex_unlock(threadArgs->mutex_);

    clock_gettime(CLOCK_REALTIME, &end);
    threadArgs->timeElapsed_ = end.tv_sec - begin.tv_sec;
    threadArgs->timeElapsed_ += (end.tv_nsec - begin.tv_nsec) / 1e9;

    pthread_exit(NULL);
}

int main(int argc, char const * const argv[])
{
    double epsilon = 0.0;
    if (2 != argc) {
        printf("Usage: %s epsilon\n", argv[0]);
        return EXIT_FAILURE;
    } else if (0 > sscanf(argv[1], "%lf", &epsilon)) {
        printf("Couldn't parse epsilon\n");
        return EXIT_FAILURE;
    }

    assert(epsilon > 0.0);

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &begin);

    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);

    thread_info_t *threads = (thread_args_t *)malloc(sizeof(thread_info_t) * N_THREADS);

    double invHighb = 1.0 / LOWB;
    double invLowb = 1.0 / HIGHB;
    double invInterval = invHighb - invLowb;

    int nSteps = ceil(sqrt(pow(invInterval, 3) / (12 * epsilon)));

    double invIntervalPerThread = invInterval / N_THREADS;
    int nStepsPerThread = nSteps / N_THREADS;

    double bank = 0;
    for (int i = 0; i < N_THREADS; ++i) {
        threads[i].args_.invBegin_ = invLowb + invIntervalPerThread * i;
        threads[i].args_.invEnd_ = invLowb + invIntervalPerThread * (i + 1);
        threads[i].args_.nSteps_ = nStepsPerThread;
        threads[i].args_.bank_ = &bank;
        threads[i].args_.mutex_ = &mutex;
    }

    for (int i = 0; i < N_THREADS; ++i) {
        pthread_create(&(threads[i].pid_), NULL, compute, &(threads[i].args_));
    }

    for (int i = 0; i < N_THREADS; ++i) {
        pthread_join(threads[i].pid_, NULL);
    }

    pthread_mutex_destroy(&mutex);

    clock_gettime(CLOCK_REALTIME, &end);
    elapsed = end.tv_sec - begin.tv_sec;
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1e9;

    printf("Final result: %lf\n", bank);

    printf("Times:\n");
    for (int i = 0; i < N_THREADS; ++i) {
        printf("\tThread [%d]: %lfs\n", i, threads[i].args_.timeElapsed_);
    }
    printf("\tMain thread: %lfs\n", elapsed);

    printf("Result:\n\tComputed: %.6lf\n\tReal: %.6lf\n\tError: %.6lf\n", bank, TRUE_VALUE, fabs(bank - TRUE_VALUE));

    free(threads);

    return EXIT_SUCCESS;
}
