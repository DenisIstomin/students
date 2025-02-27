# Черновик
# Введение
Вся теория:
* 2022 год: https://compscicenter.ru/courses/hp-course/2022-spring/classes/8861/
* 2016 год: https://compscicenter.ru/courses/hp-course/2016-spring/classes/1719/

Лабораторные работы (8 штук) можно разделить на 2 блока:
1. Низкоуровневые: программирование на языке C (SSE, Pthreads, OpenMP)
2. Высокоуровневые: программирование на языке Java (Java threads и коллекции,  Java util concurrent, Java IPC socket, Java IPC mappedbus)
Последняя лаб. работа посвященная MPI является обзорной и дополнительной.

Доп. материалы и презентации есть в отдельном архиве: todo ссылка на gdrive

# Настройка окружения

Настройка окружения:
https://code.visualstudio.com/docs/cpp/config-wsl

# Примеры кода
В каждой из лабораторных работ приведен пример кода, комментария вида
~~~c
// ToDo: что-то
~~~
означают, что некоторый код отсутствует и его необходимо написать.
# Лабораторные работы
## 1. SSE
Необходимо:
1. При помощи SSE инструкций написать программу (или функцию), которая перемножает массив из 4х чисел размером 32 бита.
2. Написать аналогичную программу (или функцию) которая решает ту же задачу последовательно.
3. Сравнить производительность
4. Проанализировать сгенерированный ассемблер: ```gcc -S sse.c```

https://felix.abecassis.me/2011/09/cpp-getting-started-with-sse/

http://www.cs.fsu.edu/~baker/opsys/notes/assembly.html

~~~c
void sse(float a[], float b[], float c[]) {
  asm volatile (
                "movups %[a], %%xmm0\n"
                "movups %[b], %%xmm1\n"
                "mulps %%xmm1, %%xmm0\n"
                "movups %%xmm0, %[c]\n"
                :
                : [a]"m"(*a), [b]"m"(*b), [c]"m"(*c)
                : "%xmm0", "%xmm1");
  for (int i = 0; i < 4; i++) {
    printf("%f ", c[i]);
  }
  printf("\n");
}

int main(int argc, char** argv) {
  // ToDo: инициализация массивов `a` и `b`

  for (int i = 0; i < iterations_num; i++) {
    sse(a, b, c);
  }

  return 0;
}
~~~

## 2. Pthreads
Необходимо:
1. При помощи Pthreads написать программу (или функцию), которая создает `n` потоков и каждый из потоков выполняет длительную операцию.
2. Написать аналогичную программу (или функцию) которая решает ту же задачу последовательно.
3. Сравнить производительность

~~~c
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

int counter = 0;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


void *heavy_task(void *i) {
  int thread_num = *((int*) i);

  // Пример критической секции с мьютексом
  printf("\tThread #%d started\n", thread_num);
  pthread_mutex_lock(&mutex);
  printf("\t\tThread #%d acquired mutex\n", thread_num);
  counter++;
  printf("\t\t\tThread #%d, counter: %d\n", thread_num, counter);
  printf("\t\tThread #%d released mutex\n", thread_num);
  pthread_mutex_unlock(&mutex);

  // Long-running task
  for (int i = 0; i < 1e8; i++) {
    sqrt(i);
  }

  printf("\tThread #%d finished\n", thread_num);
  free(i);
}

void pthreads(int threads_num) {

  pthread_t threads[threads_num];
  int status;

  for (int i = 0; i < threads_num; i++) {

    printf("MAIN: starting thread %d\n", i);

    int *thread_num = (int*) malloc(sizeof(int));
    *thread_num = i;

    status = pthread_create(&threads[i], NULL, heavy_task, thread_num);

    if (status != 0) {
      fprintf(stderr, "pthread_create failed, error code %d\n", status);
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < threads_num; i++) {
    pthread_join(threads[i], NULL);
  }
}

int main(int argc, char** argv) {
  int threads_num = atoi(argv[1]);
  pthreads(threads_num);
  pthread_mutex_destroy(&mutex);
  return 0;
}

~~~
## 3. OpenMP (сумма элементов массива, сравнение производительности с Pthreads, последовательным кодом)
Необходимо:
1. При помощи OpenMP написать программу (или функцию), которая создает `n` потоков и каждый из потоков выполняет длительную операцию.
2. Сравнить с последовательной программой и программой с Pthreads из предыдущей лабораторной работы.

~~~c
void *heavy_task() {
  int limit = 1e8;
  for (int i = 0; i < limit; i++) {
    sqrt(i);
  }
}

void openmp(int thread_num) {
  // omp_set_dynamic(0);
  // omp_set_num_threads(thread_num);
  // printf("OpenMP threads: %d\n", omp_get_num_threads());
  #pragma omp parallel for num_threads(thread_num)
  for (int i = 0; i < thread_num; i++) {
    heavy_task();
  }
}

// ToDo: дописать функцию main()
~~~

## 4. Java Threads (тестирование работы HashMap, Hashtable, synchronized HashMap, ConcurrentHashMap в многопоточном режиме)
## 5. Java util concurrent (создание собственного считающего семафора)
## 6. Java IPC (Обмен сообщениями через socket)
## 7. Java IPC (Использование библиотеки MappedBusWriter: запустить example)
## 8*. MPI (посчитать число pi)
Обзорная работа, нужно изучить код:
~~~c
/**
 * (C) 2001 by Argonne National Laboratory
 * https://github.com/pmodels/mpich/blob/master/examples/cpi.c
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>

double f(double);

double f(double a)
{
    return (4.0 / (1.0 + a * a));
}

int main(int argc, char *argv[])
{
    int n, myid, numprocs, i;
    double PI25DT = 3.141592653589793238462643;
    double mypi, pi, h, sum, x;
    double startwtime = 0.0, endwtime;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    fprintf(stdout, "Process %d of %d is on %s\n", myid, numprocs, processor_name);
    fflush(stdout);

    n = 10000;  /* default # of rectangles */
    if (myid == 0)
        startwtime = MPI_Wtime();

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    h = 1.0 / (double) n;
    sum = 0.0;
    /* A slightly better approach starts from large i and works back */
    for (i = myid + 1; i <= n; i += numprocs) {
        x = h * ((double) i - 0.5);
        sum += f(x);
    }
    mypi = h * sum;

    MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        endwtime = MPI_Wtime();
        printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
        printf("wall clock time = %f\n", endwtime - startwtime);
        fflush(stdout);
    }

    MPI_Finalize();
    return 0;
}

~~~~
