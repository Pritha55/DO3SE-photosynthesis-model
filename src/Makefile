SRCS := \
	./modules/constants/DO3SE_ModelConstants_ml.F90 \
	./modules/constants/DO3SE_PhysicalConstants_ml.F90 \
	./modules/config/DO3SE_ConfigTypes_ml.F90 \
	./modules/util/DO3SE_Util_ml.F90 \
	./modules/photosynthesis/DO3SE_Photosynthesis_mocks.F90 \
	./modules/photosynthesis/DO3SE_Photosynthesis_ml.F90
OBJS := $(SRCS:%.F90=%.o)

FFLAGS += -I./modules/constants
FFLAGS += -I./modules/config
FFLAGS += -I./modules/util
FFLAGS += -I./modules/photosynthesis

ifeq ($(TEST_MODE), 1)
all: build_modules

else
all: build_modules build
endif


build_modules:
	@echo ------------BUILDING MODULES---------------------------------------------
	$(MAKE) -C modules all

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

build: build_modules
	@echo ------------BUILDING PROGRAM---------------------------------------------
	$(FC) -c $(FFLAGS) main.f90
	# $(FC) $(FFLAGS) -o run main.f90 ./modules/module_demo/*.o ./modules/module_ext/*.o
	$(FC) $(FFLAGS) -o run main.f90 $(OBJS)
clean:
	@echo -----------------CLEANING----------------------------------------
	$(MAKE) -C modules clean
	$(RM) *.o *.mod *.a
	$(RM) run