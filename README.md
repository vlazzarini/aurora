Aurora
========

A minimal header-only C++ audio synthesis and processing toolkit.


Getting Started
------

Aurora is a collection of header files which can be included in your
C++ source code to implement audio signal processing. It provides
generic forms of the most common synthesis and transformation
components. The underlying principles of Aurora are:

- Signals are passed to and from processing objects as `std::vector`
objects holding a vector size (`vsize`) of  audio samples.

- Inputs are always read-only.

- Processing objects hold their own output signals (as `std::vector`
objects).

- Processing objects always produce a single channel of audio data.
Multichannel audio data needs to be managed by keeping separate
audio streams.

- The audio processing interface employs *functional* operators with
various types of parameters, which depend on the kind of processing
object used.  These functions always return a read-only reference to
the audio data object held by the processing object.

- The sample type is a template argument for all libraries. These are
normally expected to be a floating-number type (e.g. `float` or `double`),
although integral types may work with some objects in some cases, but
are not recommended.

- Audio-data consuming objects will always adjust their data vector
size to match their input sizes. If two or more inputs are used, the
vector size is adjusted to the length of the shortest (if lengths are
different). Fence-post errors are thus avoided, but programmers should
always ensure that vector sizes match.

- Vector sizes of objects can be changed during processing without
reallocation provided that they have enough capacity.
If necessary, this can be ensured by reserving memory by means of the `prealloc()`
method.

- Some objects provide a `reset()` method as part of their interface.
These methods should be invoked whenever the sampling rate changes.

- As with the standard C++ library, users are not supposed to derive
classes from the ones provided by Aurora, but use them to create
different types of signal processing graphs by composition.

Generic Interface
----------

Some of the classes in Aurora provide a generic interface that makes
use of processing functions supplied to objects. For example, the
`Osc` class can be used to construct various types of oscillators,
depending on the oscillator function supplied. Functions for the
most common forms of processing are supplied.
Some examples of these can be found in the
relevant header files.

Examples
-----

A variety of basic examples are provided. Some of these have no
dependencies and output ASCII floating-point samples to the standard
output. Others make use of libsndfile to access soundfiles of various
formats.

A simple usage example is given by the following code employing
a wavetable oscillator and an ADSR envelope, which are composed
into a synthesis class:

```
using namespace Aurora;
struct Synth {
  float att, dec, sus; // envelope parameters
  std::vector<float> wave; // wavetable
  Env<float> env;  // envelope object
  Osc<float, lookupi> osc;  // oscillator object

 // Synth constructor
 // rt: release time
 // sr: sampling rate
  Synth(float rt, float sr = def_sr)
      : att(0.f), dec(0.f), sus(0.f), wave(def_vsize),
        env(ads_gen(att, dec, sus), rt, sr), osc(&wave, sr) {
    std::size_t n = 0;
    for (auto &s : wave) {
      s = sin<float>((1. / wave.size()) * n++);
    }
  };

  // synth function
  // a: amplitude
  // f: frequency
  // gate: envelope gate 
  const std::vector<float> &operator()(float a, float f, bool gate) {
    return env(osc(a, f), gate);
  }
};
```

