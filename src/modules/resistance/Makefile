SRCS := DO3SE_Resistance_ml.F90
OBJS := $(SRCS:%.F90=%.o)
FFLAGS += -I../constants
FFLAGS += -I../util

all: libsut.a

libsut.a: $(OBJS)
	$(AR) -r $@ $?

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a