if [ ! -d "build" ]; then
   mkdir build;
fi

cython -a css_cython_parallel.pyx

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c css_cython_parallel.c -o build/css_cython_parallel.o

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c css.c -o build/css.o

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c comparative.c -o build/comparative.o

icc -pthread -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c threadcss.c -o build/threadcss.o

icc -pthread -shared -O3 -xAVX -mavx build/css_cython_parallel.o build/css.o build/comparative.o build/threadcss.o -L/cluster/software/VERSIONS/python2-2.7.3/lib -L//cluster/software/VERSIONS/gsl-1.16/lib -lm -lgsl -lgslcblas -lpython2.7 -o /hyperbrowser/src/hb_core_comparative/quick/statistic/css_cython_parallel.so