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

# Лабораторные работы
## 1. SSE
(вычисление квадратного корня или суммы массивов, анализ сгенерированного ассемблера)

https://felix.abecassis.me/2011/09/cpp-getting-started-with-sse/

http://www.cs.fsu.edu/~baker/opsys/notes/assembly.html

Сравнение с последовательным кодом

## 2. Pthreads (сумма элементов массива, сравнение производительности с последовательным кодом)
## 3. OpenMP (сумма элементов массива, сравнение производительности с Pthreads, последовательным кодом)
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
