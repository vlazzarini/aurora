// SndBase.h
// Base class
// Binary operations, utilities
//
// (c) V Lazzarini, 2021
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE

#ifndef _AURORA_SNDBASE_
#define _AURORA_SNDBASE_

#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace Aurora {
const int def_vsize = 64;
const double def_sr = 44100.;

/** SndBase class \n
    Aurora Library base class \n
    S: sample type
*/
template <typename S> class SndBase {
  std::vector<S> sig;

protected:
  const std::vector<S> &process(std::function<S()> f, std::size_t sz) {
    if (sz)
      vsize(sz);
    for (auto &s : sig)
      s = f();
    return sig;
  }

  std::vector<S> &get_sig() { return sig; };

public:
  /** Constructor \n
      vsize: signal vector size
  */
  SndBase(std::size_t vsize = def_vsize) : sig(vsize){};

  /** Vector size query \n
      returns the object vector size
  */
  std::size_t vsize() const { return sig.size(); }

  /** Vector size setting \n
      n: new vector size
  */
  void vsize(std::size_t n) { sig.resize(n); }

  /** Vector access \n
      returns the object vector
  */
  const std::vector<S> &vector() const { return sig; }

  /** Preallocate vector memory \n
      size: size of vector to reserve in memory \n
      this method does not change the vector size.
  */
  void prealloc(std::size_t size) { sig.reserve(size); }

  /** Copy sample data to out \n
      out: buffer receiving the audio data \n
      output buffer needs to have space for vsize() samples.
  */
  void copy_out(S *out) { std::copy(sig.begin(), sig.end(), out); }
};

/** Buff class \n
    sample buffer \n
    S: sample type
*/
template <typename S> class Buff : public SndBase<S> {
  using SndBase<S>::get_sig;

public:
  /** Constructor \n
      vsize: signal vector size
  */
  Buff(std::size_t vsize = def_vsize) : SndBase<S>(vsize){};

  /** Buffer \n
      in: audio input array \n
      n: number of samples in audio input array \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const S *in, std::size_t n = def_vsize) {
    auto &s = get_sig();
    this->resize(n);
    std::copy(in, in + n, s.begin());
    return s;
  }
};

/** BinOp class \n
    Binary operations \n
    OP: binary operation function \n
    S: sample type
*/
template <typename S, S (*OP)(S, S)> class BinOp : public SndBase<S> {
  using SndBase<S>::process;

public:
  /** Constructor \n
      vsize: signal vector size
  */
  BinOp(std::size_t vsize = def_vsize) : SndBase<S>(vsize){};

  /** Binary operation \n
      a: scalar input \n
      s: signal input \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(S a, const std::vector<S> &s) {
    std::size_t n = 0;
    return process([&]() { return OP(a, s[n++]); }, s.size());
  }

  /** Binary operation \n
      s: signal input \n
      a: scalar input \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &s, S a) {
    std::size_t n = 0;
    return process([&]() -> S { return OP(a, s[n++]); }, s.size());
  }

  /** Binary operation \n
      s1: signal input 1 \n
      s2: signal input 2 \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &s1,
                                   const std::vector<S> &s2) {
    std::size_t n = 0;
    return process(
        [&]() -> S {
          auto s = OP(s1[n], s2[n]);
          n++;
          return s;
        },
        s1.size() < s2.size() ? s1.size() : s2.size());
  }
};

/** Mix class \n
    n-signal mixer \n
    S: sample type
*/
template <typename S> class Mix : public SndBase<S> {
  using SndBase<S>::process;
  using SndBase<S>::get_sig;

  const std::vector<S> &mix(std::vector<S> &out, const std::vector<S> &in) {
    std::size_t n = 0;
    return process(
        [&]() {
          auto s = out[n] + in[n];
          n++;
          return s;
        },
        in.size() < this->vsize() ? in.size() : this->vsize());
  }

  template <typename... Ts>
  const std::vector<S> &mix(std::vector<S> &out, const std::vector<S> &in,
                            Ts... args) {
    mix(out, in);
    return mix(out, args...);
  }

public:
  /** Constructor \n
      vsize: signal vector size
  */
  Mix(std::size_t vsize = def_vsize) : SndBase<S>(vsize){};

  /** Mixer \n
      in: first input signal \n
      args: any number of other signal inputs
  */
  template <typename... Ts>
  const std::vector<S> &operator()(const std::vector<S> &in, Ts... args) {
    auto &s = get_sig();
    s.resize(in.size());
    std::fill(s.begin(), s.end(), 0);
    return mix(s, in, args...);
  }
};

/** linear interpolation table lookup \n
    S: sample type \n
    pos: reading position (no bounds check) \n
    t: table
*/
template <typename S>
inline S linear_interp(double pos, const std::vector<S> &t) {
  size_t posi = (size_t)pos;
  double frac = pos - posi;
  return t[posi] +
         frac * ((posi != t.size() - 1 ? t[posi + 1] : t[0]) - t[posi]);
}

/** cubic interpolation table lookup \n
    S: sample type \n
    pos: reading position (no bounds check) \n
    t: table
*/
template <typename S>
inline S cubic_interp(double pos, const std::vector<S> &t) {
  size_t posi = (size_t)pos;
  double frac = pos - posi;
  double a = posi == 0 ? t[0] : t[posi - 1];
  double b = t[posi];
  double c = posi != t.size() - 1 ? t[posi + 1] : t[0];
  double d =
      posi != t.size() - 2 ? (posi != t.size() - 1 ? t[posi + 1] : t[0]) : t[1];
  double tmp = d + 3.f * b;
  double fracsq = frac * frac;
  double fracb = frac * fracsq;
  return fracb * (-a - 3 * c + tmp) / 6 + fracsq * ((a + c) / 2 - b) +
         frac * (c + (-2 * a - tmp) / 6) + b;
}

} // namespace Aurora

#endif // _AURORA_SNDBASE_
