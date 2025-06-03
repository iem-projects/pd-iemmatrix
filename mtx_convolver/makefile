all: array convolver test

array: src/array.c 
	gcc -c src/array.c -o bin/array.o 

convolver: bin/array.o src/convolver.c
	gcc -c src/convolver.c -o bin/convolver.o 

test: bin/array.o bin/convolver.o src/test_convolver.c
	gcc bin/array.o bin/convolver.o src/test_convolver.c bin/libfftw3f-3.dll -o bin/test.exe

clean:
	rm bin/*.o bin/*.exe
