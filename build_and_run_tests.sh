#!/bin/bash

cd src && \
make clean && \
make all TEST_MODE=1 && \
echo "----------------------- RUNNING TESTS -----------------------------"
./modules/photosynthesis/run_tests || true &&\
make clean