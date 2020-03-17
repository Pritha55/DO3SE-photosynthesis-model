#!/bin/bash

cd src/modules && \
make clean && \
make all ./run_tests && \
./run_tests || true &&\
make clean 