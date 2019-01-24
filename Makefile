#
CC_FLAGS = -Wall   

LD_FLAGS =  -lfftw3 -lm

OBJ = cdcProt_input.o cdcProt_attach.o cdcProt_calculate.o cdcProt_main.o\

SRC = cdcProt_input.c cdcProt_attach.c cdcProt_calculate.c cdcProt_main.c\

HDR = cdcProt_head.h
 

CC = gcc     # Compiler

EXECUTABLE = cdcProt



$(EXECUTABLE) : $(SRC) $(HDR)                 #$(OBJ)
#
	$(CC) $(CC_FLAGS) -o cdcProt $(SRC) $(LD_FLAGS)
