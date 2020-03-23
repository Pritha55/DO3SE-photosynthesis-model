SRCS := DO3SE_Met_ml.F90
OBJS := $(SRCS:%.F90=%.o)

FFLAGS += -I../constants
FFLAGS += -I../config
FFLAGS += -I../resistance
FFLAGS += -I../util

all: libsut.a

libsut.a: $(OBJS)
	$(AR) -r $@ $?

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a