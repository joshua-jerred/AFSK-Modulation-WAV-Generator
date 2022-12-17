CC = g++
CFLAGS = #-Wall -Wpedantic -Wextra -Werror

#all: main
#
#main: main.o afsk.o ax25.o
#	$(CC) $(CFLAGS) -o main main.o afsk.o
#
#ax25.o: ax25.cpp ax25.h
#	$(CC) $(CFLAGS) -c ax25.cpp

#afsk.o: afsk.cpp afsk.h
#	$(CC) $(CFLAGS) -c afsk.cpp

afsk: afsk.cpp afsk.h
	$(CC) $(CFLAGS) -o afsk afsk.cpp

clean:
	rm -f *.o *.wav