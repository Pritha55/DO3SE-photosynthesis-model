SRCS := \
	../constants/DO3SE_ModelConstants_ml.F90 \
	../constants/DO3SE_PhysicalConstants_ml.F90 \
	../config/DO3SE_ConfigTypes_ml.F90 \
	../util/DO3SE_Util_ml.F90 \
	../resistance/DO3SE_Resistance_ml.F90 \
	DO3SE_Photosynthesis_mocks.F90 \
	DO3SE_Photosynthesis_helpers_ml.F90 \
	DO3SE_Photosynthesis_far_ml.F90 \
	DO3SE_Photosynthesis_ml.F90

OBJS := $(SRCS:%.F90=%.o)



all: libsut.a


# Add any compile flags?
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I../photosynthesis
FFLAGS += -I../constants
FFLAGS += -I../config
FFLAGS += -I../util
FFLAGS += -I../resistance

# Setup/compile library?
libsut.a: $(OBJS)
	$(AR) -r $@ $?


# Tells make to Compile using FC compiler
%.o : %.F90
	$(FC) -c $(FFLAGS) $<

# ==== ======================================================DEFINE TESTS ==== #
# Define Tests and their dependencies
tests: photosynthesis_tests
# Get pFUnit
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

OTHER_LIBRARIES = -L../photosynthesis -L../constants -L../config -L../util -L../resistance -lsut

# **run_tests**
photosynthesis_tests_TESTS := DO3SE_Photosynthesis_test.pf DO3SE_Photosynthesis_far_test.pf DO3SE_Photosynthesis_helpers_test.pf
photosynthesis_tests_OTHER_LIBRARIES := $(OTHER_LIBRARIES)

$(eval $(call make_pfunit_test,photosynthesis_tests))


# ==== CLEAN UP ==== #
clean:
	$(RM) *.o *.mod *.a *.inc *_tests *_test.F90
