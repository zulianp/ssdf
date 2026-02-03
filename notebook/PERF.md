# Performance of Unsigned Distance Function

Results with CUBIQL integration with support for GPU offloading

## Bone geometry

13029862 query   points 
11684912 surface elements

NVIDIA Hopper GPU

1  GPU  4.17183  [s]
16 GPUs 0.701588 [s] (6x speed-up, 37% strong scaling efficiency)

## Pump geometry

316200 query   points
110244 surface elements

Arm M1 8 cores

2.36267 [s]