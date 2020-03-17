SRCS := DO3SE_Util_ml.F90
OBJS := $(SRCS:%.F90=%.o)
FFLAGS += -I../constants

all: libsut.a

libsut.a: $(OBJS)
	$(AR) -r $@ $?

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a