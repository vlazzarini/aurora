// Osc.h:
// Generic oscillator
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

#ifndef _AURORA_OSC_
#define _AURORA_OSC_

#include "SndBase.h"
#include <cmath>
#include <functional>

namespace Aurora {
const double twopi = 2 * M_PI;

/** Table lookup function generator for Osc \n
    S: sample type \n
    t: function table  \n
    returns a truncating table lookup function
*/
template <typename S> std::function<S(S)> lookup_gen(const std::vector<S> &t) {
  return [&t](double ph) -> S { return t[(std::size_t)(ph * t.size())]; };
}

/** Table lookup function generator for Osc \n
    S: sample type \n
    t: function table  \n
    returns an interpolating table lookup function
*/
template <typename S> std::function<S(S)> lookupi_gen(const std::vector<S> &t) {
  return [&t](double ph) -> S { return linear_interp(ph * t.size(), t); };
}

/** Sine function for Osc \n
    S: sample type \n
    ph: normalised phase  \n
    returns the sine of ph*2*$M_PI
*/
template <typename S> S sin(double ph) { return (S)std::sin(ph * twopi); }

/** Cosine function for Osc \n
    S: sample type \n
    ph: normalised phase \n
    returns the cosine of ph*2*$M_PI
*/
template <typename S> S cos(double ph) { return (S)std::cos(ph * twopi); }

/** Phase function for Osc \n
    S: sample type \n
    ph: normalised phase \n
    returns ph
*/
template <typename S> S phase(double ph) { return (S)ph; }

/** Osc class  \n
    Generic oscillator \n
    S: sample type
*/
template <typename S> class Osc : public SndBase<S> {
  using SndBase<S>::process;

protected:
  double ph;
  S ts;
  std::function<S(S)> fun;

  virtual S synth(S a, S f, double &phs) {
    S s = (S)(a * fun(phs));
    phs += f * ts;
    while (phs < 0)
      phs += 1.;
    while (phs >= 1.)
      phs -= 1.;
    return s;
  }

public:
  /** Constructor \n
      f: oscillator function \n
      fs: sampling rate \n
      vsize: signal vector size
  */
  Osc(std::function<S(S)> f = cos<S>, S fs = (S)def_sr,
      std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), ph(0.), ts(1 / fs), fun(f){};

  virtual ~Osc(){};

  /** Sampling rate query \n
      returns sampling rate
  */
  S fs() const { return 1 / ts; }

  /** Oscillator \n
      a: scalar amplitude \n
      f: scalar frequency \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(S a, S f) {
    double phs = ph;
    const std::vector<S> &s = process([&]() { return synth(a, f, phs); }, 0);
    ph = phs;
    return s;
  }

  /** Oscillator \n
      a: scalar amplitude \n
      fm: frequency signal \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(S a, const std::vector<S> &fm) {
    double phs = ph;
    std::size_t n = 0;
    const std::vector<S> &s =
        process([&]() { return synth(a, fm[n++], phs); }, 0);
    ph = phs;
    return s;
  }

  /** Oscillator \n
      am: amplitude signal \n
      f: scalar frequency \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &am, S f) {
    double phs = ph;
    std::size_t n = 0;
    const std::vector<S> &s =
        process([&]() { return synth(am[n++], f, phs); }, 0);
    ph = phs;
    return s;
  }

  /** Oscillator \n
      am: amplitude signal \n
      fm: frequency signal \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &am,
                                   const std::vector<S> &fm) {
    double phs = ph;
    std::size_t n = 0;
    const std::vector<S> &s = process(
        [&]() {
          auto s = synth(am[n], fm[n], phs);
          n++;
          return s;
        },
        am.size() < fm.size() ? am.size() : fm.size());
    ph = phs;
    return s;
  }

  /** set the oscillator function
      f: oscillator function to be used
  */
  void func(std::function<S(S)> f) { fun = f; }
};
} // namespace Aurora

#endif // _AURORA_OSC_
