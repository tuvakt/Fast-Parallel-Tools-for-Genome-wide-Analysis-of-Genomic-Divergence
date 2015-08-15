if [ ! -d "build" ]; then
   mkdir build;
fi

cython -a fisher_cython_parallel.pyx

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c fisher_cython_parallel.c -o build/fisher_cython_parallel.o

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c cFisher.c -o build/cFisher.o

icc -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c comparative.c -o build/comparative.o

icc -pthread -DNDEBUG -g -fwrapv -03 -Wall -Wstrict-prototypes -xAVX -mavx -fPIC -I/cluster/software/VERSIONS/python2-2.7.3/include/python2.7 -c threadfisher.c -o build/threadfisher.o

icc -pthread -shared -O3 -xAVX -mavx build/fisher_cython_parallel.o build/cFisher.o build/comparative.o build/threadfisher.o -L/cluster/software/VERSIONS/python2-2.7.3/lib -lm -lpython2.7 -o /hyperbrowser/src/hb_core_comparative/quick/statistic/fisher_cython_parallel.so