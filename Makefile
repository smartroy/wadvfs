cc=g++
cc2=gcc
test_core:Core.o test.o flp.o util.o npe.o shape.o RCutil.o tilts.o temperature.o
	$(cc) -o  test_core Core.o test.o flp.o util.o npe.o shape.o RCutil.o tilts.o temperature.o temperature_block.o temperature_grid.o
test.o:test.cc Core.h Core.o
	$(cc) -c test.cc
Core.o:Core.h Core.cc flp.o
	$(cc) -c Core.cc -std=c++11
temperature.o:temperature.c temperature.h temperature_block.o temperature_grid.o temperature_block.h temperature_grid.h 
	$(cc2)  -DDEBUG -O0 -g -Wall -c temperature.c -ltmperature_block -ltemperature_grid
temperature_block.o:temperature_block.c temperature_block.h
	$(cc2) -DDEBUG -O0 -g -Wall -c temperature_block.c  -lm
temperature_grid.o:temperature_grid.c temperature_grid.h
	$(cc2) -DDEBUG -O0 -g -Wall -c temperature_grid.c -lm
tilts.o:tilts.c tilts.h
	$(cc2)  -c -DDEBUG -O0 -g -Wall tilts.c -lm
RCutil.o:RCutil.c util.o
	$(cc2) -DDEBUG -O0 -g -Wall -c RCutil.c -lm
flp.o:flp.c flp.h util.h npe.o shape.o
	$(cc2)  -DDEBUG -O0 -g -Wall -c flp.c -lm -lnpe -lshape
util.o:util.c util.h
	$(cc2)  -DDEBUG -O0 -g -Wall -c util.c -lm
npe.o:npe.c npe.h
	$(cc2) -DDEBUG -O0 -g -Wall -c npe.c -lm
shape.o:shape.c shape.h
	$(cc2) -DDEBUG -O0 -g -Wall -c shape.c -lm

clean:
	rm *.o
