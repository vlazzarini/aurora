// Env.h:
// Generic envelope
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

#ifndef _AURORA_ENV_
#define _AURORA_ENV_

#include "SndBase.h"
#include <functional>

namespace Aurora {

/** ADS function generator for Env \n
  S: sample type \n
  a: attack  \n
  d: decay  \n
  s: sustain \n
  returns an ADS envelope function
*/
template <typename S>
std::function<S(double, S, S)> ads_gen(const S &a, const S &d, const S &s) {
  return [&a, &d, &s](double t, S e, S ts) -> S {
    if (t < a && e < 1)
      return e + ts / a;
    else if (t < a + d)
      return e + (s - 1) * ts / d;
    else
      return s;
  };
}

/** Env class  \n
    Generic envelope \n
    S: sample type
*/
template <typename S> class Env : public SndBase<S> {
  using SndBase<S>::process;
  std::function<S(double, S, S)> fun;
  double time;
  S prev;
  S ts;
  S fac;

  S synth(S e, double &t, bool gate, std::function<S(double, S, S)> fn) {
    S s;
    if (gate) {
      s = fn(t, e, ts);
      t += ts;
    } else {
      if (e < 0.00001)
        s = 0;
      else
        s = e * fac;
      t = 0;
    }
    return s;
  }

public:
  /** Constructor \n
      f:  envelope function \n
      rt: release time \n
      fs: sampling rate \n
      vsize: signal vector size
  */
  Env(std::function<S(double, S, S)> f, S rt, S fs = def_sr,
      std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), fun(f), prev(0), ts(1. / fs),
        fac(std::pow(0.001, ts / rt)){};

  /** Release time setter
     rt: release time
   */
  void release(S rt) { fac = std::pow(0.001, ts / rt); }

  /** Envelope \n
      gate: envelope gate
   */
  const std::vector<S> &operator()(bool gate) {
    double t = time;
    S e = prev;
    auto &s = process([&]() { return (e = synth(e, t, gate, fun)); }, 0);
    prev = e;
    time = t;
    return s;
  }

  /** Envelope \n
      offs: sig offset
      scal: sig scale
      gate: envelope gate
   */
  const std::vector<S> &operator()(S offs, S scal, bool gate) {
    double t = time;
    S e = prev;
    auto &s = process(
        [&]() {
          e = synth(e, t, gate, fun);
          return (e * scal + offs);
        },
        0);
    prev = e;
    time = t;
    return s;
  }

  /** Envelope \n
      a: input signal
      gate: envelope gate
   */
  const std::vector<S> &operator()(const std::vector<S> &a, bool gate) {
    double t = time;
    S e = prev;
    std::size_t n = 0;
    auto &s = process(
        [&]() {
          e = synth(e, t, gate, fun);
          return e * a[n++];
        },
        a.size());
    prev = e;
    time = t;
    return s;
  }

  /** Sampling rate query \n
      returns sampling rate
  */
  S fs() const { return 1 / ts; }

  /** reset the envelope \n
     fs: sampling rate
  */
  void reset(S fs) {
    time = 0;
    prev = 0;
    ts = 1 / fs;
  }

  /** set the envelope function \n
      f: envelope function to be used
   */
  void func(const std::function<S(double, S, S)> f) { fun = f; }
};
} // namespace Aurora

#endif // _AURORA_ENV_
