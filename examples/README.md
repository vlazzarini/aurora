Aurora Examples
=================

Examples can be built using

```
c++ -o example example.cpp -I../include -std=c++14
```

**wave.cpp**: classic bandlimited waveform  
**custom.cpp**: custom bandlimited waveform  
**stackedfm.cpp**: stacked frequency modulation

In addition, there is a utility to convert from ASCII samples to RIFF-Wave files. This
requires libsndfile to be installed:

```
cc -o towav towav.c -lsndfile
```

With this, it is possible to produce a soundfile by piping the standard
output of the example programs, e.g.

```
./wave 1 1 120 | ./towav out.wav
```

