  # the compiler: gcc for C program, define as g++ for C++
  CC = g++
  CFLAGS  = -g -Wall -std=c++11
   EXTRAFLAGS = -fopenmp -O3 -march=native -mfpmath=sse
 # EXTRAFLAGS = 
   
  # the build target executable:
  TARGET = main
  SOURCES = main.cpp basis_solver.cpp
  
  all: $(TARGET)

  $(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(EXTRAFLAGS) -o $(TARGET) $(SOURCES)

  clean:
	$(RM) $(TARGET)
#make clean; make; valgrind --leak-check=yes --track-origins=yes ./main > vallog.txt 2>&1 