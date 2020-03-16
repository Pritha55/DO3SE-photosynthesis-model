SRCS := mocked_functions.F90 DO3SE_ModelConstants_ml.F90 DO3SE_Util_ml.F90 DO3SE_ConfigTypes_ml.F90  DO3SE_PhysicalConstants_ml.F90 DO3SE_Resistance_ml.F90 DO3SE_Photosynthesis_ml.F90
OBJS := $(SRCS:%.F90=%.o)

# Get pFUnit
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

all: libsut.a run_tests 


# Add any compile flags?
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I./

# Setup/compile library?
libsut.a: $(OBJS)
	$(AR) -r $@ $?

# Set build order - this forces module_demo_external to be built before module_demo
# DO3SE_Photosynthesis_ml.o: module_demo_external.mod
# module_demo_external.mod: module_demo_external.o


# Tells make to Compile using FC compiler
%.o : %.F90
	$(FC) -c $(FFLAGS) $<


# ==== DEFINE TESTS ==== #
# Define Tests and their dependencies

# **run_tests**
run_tests_TESTS := DO3SE_Photosynthesis_test.pf
# run_tests_REGISTRY :=
# run_tests_OTHER_SOURCES := mocked_functions.F90 DO3SE_PhysicalConstants_ml.F90
run_tests_OTHER_LIBRARIES := -L. -lsut
# run_tests_OTHER_INCS :=
# run_tests_OTHER_SRCS := mocked_functions.F90 DO3SE_PhysicalConstants_ml.F90
# run_tests_EXTRA_USE := mocked_functions_mod DO3SE_PhysicalConstants_ml.mod
 
$(eval $(call make_pfunit_test,run_tests))


# ==== CLEAN UP ==== #
clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) DO3SE_Photosynthesis_test.F90
	${RM} run_tests
