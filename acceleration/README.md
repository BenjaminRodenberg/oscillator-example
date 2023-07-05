`mass-spring.py` solves mass-spring-system in partitioned fashion.


# Acceleration case

The goal of this case is to investigate the convergence of different acceleration schemes based on preCICE v2.5.0.2. This case builds up on the multirate case.

## Config file

See description in multirate case. Additionally, one can specify the acceleration scheme being used.

## Running

```
python3 mass-spring.py Mass-Left
```

and

```
python3 mass-spring.py Mass-Right
```

## Experiments

TODO