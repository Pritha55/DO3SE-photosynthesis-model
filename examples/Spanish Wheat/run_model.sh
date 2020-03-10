#!/bin/bash

../../tools/run_model.sh \
    "config.nml" \
    "dummy, dummy, dd, hr, met%Ts_C, met%VPD, met%u, met%precip, met%P, met%O3, met%PPFD" \
    "Spanish data.csv" \
    "this%MLMC(1,1)%leaf_gsto, this%MLMC(1,1)%bulk_gsto" \
    "output.csv"
