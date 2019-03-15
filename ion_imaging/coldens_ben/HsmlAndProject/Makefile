# PERIODIC BOUNDARY CONDITIONS:
# Compile with first two comments below removed.

EXEC  = HsmlAndProject.so
#EXEC  = HsmlAndProject_perbc.so

OBJS   = main.o tree.o peano.o allvars.o kernel.o

#OPTIONS += -DPERIODIC

# KERNEL SHAPE:
# Choose one of -DSPH_KERNEL_GADGET, -DSPH_KERNEL_CUBIC,
# -DSPH_KERNEL_QUARTIC, -DSPH_KERNEL_QUINTIC, 
# -DSPH_KERNEL_C2, -DSPH_KERNEL_C4, -DSPH_KERNEL_C6

KERNEL = -DSPH_KERNEL_C2

INCL   = proto.h allvars.h kernel.h Makefile

CFLAGS =  -O3 $(OPTIONS) $(KERNEL) -fPIC -Wall -mcmodel=medium
LNKCMD =  ld -shared

LIBS   =  -lm -g

CC     =  gcc

$(EXEC): $(OBJS)
	 $(LNKCMD)  $(OBJS) $(LIBS)  -o $(EXEC)

$(OBJS): $(INCL)

clean:
	 rm -f $(OBJS) $(EXEC)
