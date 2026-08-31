#!/bin/bash

rsync -avz -e ssh ./ take@10.35.22.56:/home/take/TMTGEMv1.3/tmtgem/solver/
