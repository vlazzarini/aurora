Utility Programs
=================
These utilities require libsndfile.

towav
----

**towav.c**: conversion from ASCII decimal samples at stdin to RIFF-Wave files. 


```
cc -o towav towav.c -lsndfile
```

Usage: 

```
source |  towav out.wav [sr] [channels]
```

where `source` is a program that produces decimal samples (-1. to 1.)
to stdin.

wtgen
-------

Wavetable generation from a soundfile to C++ header file.

```
cc -o wtgen wtgen.c -lsndfile
```

Usage:

```
wtgen infile outfile namespace start(samples) end(samples) basefr
```

where `infile` is a soundfile and `outfile` is the C++ header file
name. The data is put to a `std::vector<float>` inside the `namespace`
given. The wavetable is taken from `start` and `end` in samples,
and a base frequency for the waveform should be given in `basefr`.
