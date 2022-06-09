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
const int def_ftlen = 16384;

/** Truncating table lookup for Osc \n
    S: sample type \n
    ph: phase \n
    t: function table  \n
    returns a sample
*/
template <typename S> inline S lookup(double ph, const std::vector<S> *t) {
  return (*t)[(std::size_t)(ph * t->size())];
}

/** Linear interp table lookup function for Osc \n
    S: sample type \n
    ph: phase  \mn
    t: function table  \n
    returns an interpolated sample
*/
template <typename S> inline S lookupi(double ph, const std::vector<S> *t) {
  return linear_interp(ph * t->size(), *t);
}

/** Cubic interp table lookup function for Osc \n
    S: sample type \n
    ph: phase  \mn
    t: function table  \n
    returns an interpolated sample
*/
template <typename S> S inline lookupc(double ph, const std::vector<S> *t) {
  return cubic_interp(ph * t->size(), *t);
}

/** Sine function for Osc \n
    S: sample type \n
    ph: normalised phase  \n
    returns the sine of ph*2*$M_PI
*/
template <typename S> inline S sin(double ph, const std::vector<S> *t = 0) {
  (void)t;
  return (S)std::sin(ph * twopi);
}

/** Cosine function for Osc \n
    S: sample type \n
    ph: normalised phase \n
    returns the cosine of ph*2*$M_PI
*/
template <typename S> inline S cos(double ph, const std::vector<S> *t = 0) {
  (void)t;
  return (S)std::cos(ph * twopi);
}

/** Phase function for Osc \n
    S: sample type \n
    ph: normalised phase \n
    returns ph
*/
template <typename S> inline S phase(double ph, const std::vector<S> *t = 0) {
  (void)t;
  return (S)ph;
}

/** Osc class  \n
    Generic oscillator \n
    S: sample type \n
    FN: oscillator function
*/
template <typename S, S (*FN)(double, const std::vector<S> *) = cos>
class Osc : public SndBase<S> {
  using SndBase<S>::process;

protected:
  double ph;
  S ts;
  const std::vector<S> *tab;

  virtual S synth(S a, S f, double &phs, const std::vector<S> *t, S pm = 0) {
    phs += pm;
    while (phs < 0)
      phs += 1.;
    while (phs >= 1.)
      phs -= 1.;
    S s = (S)(a * FN(phs, t));
    phs = f * ts + phs - pm;
    return s;
  }

  virtual void reset_obj(S fs) {
    ts = 1 / fs;
    ph = 0.;
  }

public:
  /** Constructor \n
      t: function table
      f: oscillator function \n
      fs: sampling rate \n
      vsize: signal vector size
  */
  Osc(const std::vector<S> *t, S fs = (S)def_sr, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), ph(0.), ts(1 / fs), tab(t){};

  /** Constructor \n
    fs: sampling rate \n
    vsize: signal vector size
*/
  Osc(S fs = (S)def_sr, std::size_t vsize = def_vsize)
      : Osc(nullptr, fs, vsize){};

  virtual ~Osc(){};

  /** Sampling rate query \n
      returns sampling rate
  */
  S fs() const { return 1 / ts; }

  /** Oscillator \n
      a: scalar amplitude \n
      f: scalar frequency \n
      pm: scalar phase \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(S a, S f, S pm = 0) {
    double phs = ph;
    auto &s = process([&]() { return synth(a, f, phs, tab, pm); }, 0);
    ph = phs;
    return s;
  }

  /** Oscillator \n
      a: scalar amplitude \n
      fm: frequency signal \n
      pm: scalar phase \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(S a, const std::vector<S> &fm, S pm = 0) {
    double phs = ph;
    std::size_t n = 0;
    auto &s =
        process([&]() { return synth(a, fm[n++], phs, tab, pm); }, fm.size());
    ph = phs;
    return s;
  }

  /** Oscillator \n
      am: amplitude signal \n
      f: scalar frequency \n
      pm: scalar phase \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &am, S f, S pm = 0) {
    double phs = ph;
    std::size_t n = 0;
    auto &s =
        process([&]() { return synth(am[n++], f, phs, tab, pm); }, am.size());
    ph = phs;
    return s;
  }

  /** Oscillator \n
      am: amplitude signal \n
      fm: frequency signal \n
      pm: scalar phase \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &am,
                                   const std::vector<S> &fm, S pm = 0) {
    double phs = ph;
    std::size_t n = 0;
    auto &s = process(
        [&]() {
          auto s = synth(am[n], fm[n], phs, tab, pm);
          n++;
          return s;
        },
        am.size() < fm.size() ? am.size() : fm.size());
    ph = phs;
    return s;
  }

  /** Oscillator \n
    a: scalar amplitude \n
    f:  scalar frequency  \n
    pm: phase modulation signal \n
    returns reference to object signal vector
 */
  const std::vector<S> &operator()(S a, S f, const std::vector<S> &pm) {
    double phs = ph;
    std::size_t n = 0;
    auto &s =
        process([&]() { return synth(a, f, phs, tab, pm[n++]); }, pm.size());
    ph = phs;
    return s;
  }

  /** Oscillator \n
      am: amplitude signal \n
      f:  scalar frequency  \n
      pm: phase modulation signal \n
      returns reference to object signal vector
  */
  const std::vector<S> &operator()(const std::vector<S> &am, S f,
                                   const std::vector<S> &pm) {
    double phs = ph;
    std::size_t n = 0;
    auto &s = process(
        [&]() {
          auto s = synth(am[n], f, phs, tab, pm[n]);
          n++;
          return s;
        },
        am.size() < pm.size() ? am.size() : pm.size());
    ph = phs;
    return s;
  }

  /** set the Osc function table \n
     t: function table
 */
  void table(const std::vector<S> *t) { tab = t; }

  /** set the internal oscillator phase \n
      phs: phase
  */
  void phase(double phs) { ph = phs; }

  /** reset the oscillator \n
      fs: sampling rate
   */
  void reset(S fs) { reset_obj(fs); }
};
} // namespace Aurora

#endif // _AURORA_OSC_
