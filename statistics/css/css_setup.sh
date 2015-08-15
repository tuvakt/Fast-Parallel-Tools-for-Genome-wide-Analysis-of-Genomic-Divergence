if [ ! -d "build" ]; then
   mkdir build;
fi

cython -a css_cython.pyx

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c css_cython.c -o build/css_cython.o

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c css.c -o build/css.o

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c comparative.c -o build/comparative.o

icc -pthread -shared -O3 -xAVX -mavx build/css_cython.o build/css.o build/comparative.o -L/cluster/software/VERSIONS/python2-2.7.3/lib -L//cluster/software/VERSIONS/gsl-1.16/lib -lm -lgsl -lgslcblas -lpython2.7 -o /hyperbrowser/src/hb_core_comparative/quick/statistic/css_cython.so