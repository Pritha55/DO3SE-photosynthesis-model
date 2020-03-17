SRCS := DO3SE_ModelConstants_ml.F90 \
		DO3SE_PhysicalConstants_ml.F90
OBJS := $(SRCS:%.F90=%.o)

all: libsut.a

libsut.a: $(OBJS)
	$(AR) -r $@ $?

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a