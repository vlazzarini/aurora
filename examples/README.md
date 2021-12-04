Aurora Examples
===

Examples can be built using CMake at the top-level sources
directory. They can also be built individually as shown below.
    

ASCII output
---

These examples write ASCII floats to stdout and can be built using

```
c++ -o example example.cpp - I../include -std=c++ 14
```

**wave.cpp** : classic bandlimited waveform

**custom.cpp** : custom bandlimited waveform

**stackedfm.cpp** : stacked frequency modulation

**oscil.cpp**: table lookup oscillator and envelope

**lpwave.cpp** : lowpass filter,bandlimited oscillator, and envelope

**drive.cpp** : nonlinear distortion

**svfdrive.cpp** : nonlinear svf

**lopass.cpp** : first-order lowpass filter

Soundfile output
---

These examples read and write soundfiles.They depend on libsndfile and can be built using


```
c++ -o example example.cpp -I../ include -std = c++ 14 - lsndfile
```

**filter.cpp** : lowpass filter example

**delay.cpp** : fixed comb filter example

**reverb.cpp**: convolution reverb

