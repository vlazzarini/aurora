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
objects holding audio samples.

- Inputs are always read-only.

- Processing objects hold their own output signals (as `std::vector`
objects).

- Processing objects always produce a single channel of audio data.
Multichannel audio data needs to be managed by keeping separate
audio streams.

- The audio processing interface employs functional operators with
various types of parameters, which depend on the kind of processing
object used.  These functions always return  a read-only reference to
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
If necessary, this can be ensured by reserving memory through the `prealloc()`
method.

- Some objects provide a `reset()` method as part of their interface.
These methods should be invoked whenever the sampling rate changes.

Generic Interface
----------

Some of the classes in Aurora provide a generic interface that makes
use of processing functions supplied to objects. For example, the
Osc class can be used to construct various types of oscillators,
depending on the oscillator function employed. Functions for the
most common forms of objects are supplied. User-defined functions
can also be employed to extend these classes, also with the use of lambda
functions and closures. Some examples of these can be found in the
relevant header files.

Examples
-----

A variety of basic examples are provided. Some of these have no
dependencies and output ASCII floating-point samples to the standard
output. Others make use of libsndfile to access soundfiles of various
formats.




