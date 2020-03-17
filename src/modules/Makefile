
LATEST_PFUNIT_DIR := $(lastword $(shell echo $(wildcard $(PFUNIT_DIR)/PFUNIT-4.*) | xargs -n1 | sort -V))
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

ifeq ($(TEST_MODE), 1)
all: build testing

else
all: build
endif

build:
	@echo -------------BUILDING MODULES--------------------------------------------
	$(MAKE) -C constants all
	$(MAKE) -C util all
	$(MAKE) -C config all
	$(MAKE) -C resistance all
	$(MAKE) -C photosynthesis all
	@echo ---------------------------------------------------------

testing: build
	@echo -------------BUILDING TESTS--------------------------------------------
	# $(MAKE) -C tests all
	$(MAKE) -C photosynthesis run_tests

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a *.inc
	$(MAKE) -C constants clean
	$(MAKE) -C util clean
	$(MAKE) -C config clean
	$(MAKE) -C resistance clean
	$(MAKE) -C photosynthesis clean
	$(MAKE) -C tests clean
