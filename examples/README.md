Aurora Examples
===

Examples can be built using CMake at the top-level sources
directory. They can also be built individually as shown below.
    

ASCII output
---

These synthesis examples write ASCII floats to stdout and can be built using

```
c++ -o example example.cpp - I../include -std=c++ 14
```

**wave.cpp** : classic bandlimited waveform

**custom.cpp** : custom bandlimited waveform

**stackedfm.cpp** : stacked frequency modulation

**stackedpm.cpp** : stacked phase modulation

**oscil.cpp**: table lookup oscillator and envelope

**lpwave.cpp**: lowpass filter,bandlimited oscillator, and envelope

**drive.cpp**: nonlinear distortion

**svfdrive.cpp**: nonlinear svf

**lopass.cpp**: first-order lowpass filter

**karplus.cpp**: string physical model

**operatorpm.cpp**: operator-based phase modulation

**pwm.cpp**: variable pulse-width square wave

**noise.cpp**: white noise genarator 

**buffer.cpp**: circular buffer test

**grsynth.cpp**: granular synthesis

ASCII floats can be converted to soundfiles using one of the utility
programs provided in the **utilities** folder.

Soundfile input/output
---

These examples read and write soundfiles.They depend on libsndfile and can be built using


```
c++ -o example example.cpp -I../ include -std = c++ 14 - lsndfile
```

**filter.cpp** : lowpass filter

**delay.cpp** : fixed comb filter

**reverb.cpp**: convolution reverb

**flanger.cpp**: flanger effect

**chorus.cpp**: chorus effect

**equaliser.cpp**: graphic equaliser

**follow.cpp**: envelope follower

**freverb.cpp**: extended Schroeder reverb

**resonator.cpp**: resonator filter

**grproc.cpp**: granular processing

**freqshift.cpp**: frequency shifting

**tvconv.cpp**: time-varying convolution



