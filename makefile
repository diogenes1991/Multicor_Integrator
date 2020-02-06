CPP 	= g++
NVD 	= nvcc

CFLAGS	= -std=c++11 -O3
LIBS	= -lpthread

NFLAGS = -I/home/diogenes1991/NVIDIA_CUDA-10.2_Samples/common/inc

SOURCE 	= Test
OUT 	= Multicore_tests


all: compile

compile: 
	$(CPP) $(CFLAGS) $(SOURCE).cpp -o $(OUT)_cpp $(LIBS)
	$(NVD) $(NFLAGS) $(SOURCE).cu  -o $(OUT)_nvd
	@echo 'Compilation Succesful'
	
clean:
	@rm -f *.o *~ *.mod *# *.str
	
run: 
	./$(OUT)_cpp > $(OUT)_cpp.txt &
	./$(OUT)_nvd > $(OUT)_nvd.txt &

open:
	@kate $(OUT)_cpp.txt &
	@kate $(OUT)_nvd.txt &
	
