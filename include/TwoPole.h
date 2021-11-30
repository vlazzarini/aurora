// TwoPole.h:
// 2-pole state variable filter with optional nonlinearity
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

#ifndef _AURORA_TWOPOLE_
#define _AURORA_TWOPOLE_
#include "SndBase.h"
#include <cmath>
#include <functional>

namespace Aurora {
enum : uint32_t { HP = 0, BP, BR, LP };

/** TwoPole class  \n
    2-pole state-variable filter \n
    S: sample type
*/
template <typename S> class TwoPole : SndBase<S> {
  using SndBase<S>::sig;
  S Y[3];
  double D[2];
  S W, Fac;
  S ff, dd;
  double piosr;
  std::function<S(S)> fun;

  S filter(S in, S *y, double *s, double w, double fac, double d, S drv) {
    S lp;
    y[HP] = (in - (d + w) * s[0] - s[1]) * fac;
    S u = w * fun(y[0] * drv) * 1 / drv;
    y[BP] = u + s[0];
    s[0] = y[BP] + u;
    u = w * fun(y[1] * drv) * 1 / drv;
    lp = u + s[1];
    s[1] = lp + u;
    y[BR] = y[HP] + lp;
    return lp;
  }

  void coeffs(S f, S d) {
    S w = std::tan(f * piosr);
    Fac = 1. / (1. + w * d + w * w);
    ff = f;
    dd = d;
    W = w;
  }

public:
  /** Constructor \n
   f: overdrive function \n
   sr: sampling rate \n
   vsize: vector size
  */
  TwoPole(std::function<S(S)> f, S fs = def_sr, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), Y{0}, D{0}, W(0), Fac(0), ff(0), dd(0),
        piosr(M_PI / fs), fun(f){};

  /** Constructor \n
   sr: sampling rate \n
   vsize: vector size
  */
  TwoPole(S fs = def_sr, std::size_t vsize = def_vsize)
      : TwoPole([](S x) -> S { return x; }, fs, vsize){};

  /** Filter \n
     in: input \n
     f: frequency \n
     d: damping factor (Q reciprocal) \n
     drv: overdrive amount \n
     o: output type (LP,HP,BP,BR)
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S d, S drv,
                                   uint32_t o = LP) {
    if (f != ff || d != dd)
      coeffs(f, d);
    S w = W, fac = Fac;
    std::size_t n = 0;
    drv += 1;
    for (auto &s : sig) {
      s = filter(in[n++], Y, D, w, fac, d, drv);
      if (o <= BR)
        s = Y[o];
    }
    return sig;
  }

  /** Filter \n
  in: input \n
  f: frequency  \n
  d: damping factor (Q reciprocal) \n
  drv: overdrive amount \n
  o: output type (LP,HP,BP,BR)
 */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S d, S drv,
                                   S m) {
    if (f != ff || d != dd)
      coeffs(f, d);
    S w = W, fac = Fac;
    std::size_t n = 0;
    drv += 1;
    for (auto &s : sig) {
      s = filter(in[n++], Y, D, w, fac, d, drv) * (1 - m);
      s += Y[HP] * m;
    }
    return sig;
  }

  /** Filter \n
   in: input \n
   f: frequency  \n
   d: damping factor (Q reciprocal) \n
   drv: overdrive amount \n
   m: lowpass - band reject - highpass mix (0-1) \n
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f, S d, S drv,
                                   uint32_t o = LP) {
    S &w = W, &fac = Fac;
    std::size_t n = 0;
    drv += 1;
    for (auto &s : sig) {
      if (f[n] != ff || d != dd)
        coeffs(f[n], d);
      s = filter(in[n++], Y, D, w, fac, d, drv);
      if (o <= BR)
        s = Y[o];
    }
    return sig;
  }

  /** Filter \n
  in: input \n
  f: frequency  \n
  d: damping factor (Q reciprocal) \n
  drv: overdrive amount \n
  m: lowpass - band reject - highpass mix (0-1) \n
 */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f, S d, S drv, S m) {
    S &w = W, &fac = Fac;
    std::size_t n = 0;
    drv += 1;
    for (auto &s : sig) {
      if (f[n] != ff || d != dd)
        coeffs(f[n], d);
      s = filter(in[n++], Y, D, w, fac, d, drv) * (1 - m);
      s += Y[HP] * m;
    }
    return sig;
  }
};
} // namespace Aurora

#endif //_AURORA_TWOPOLE_
