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
// POSSIBILITY OF SUCH DAMAGE

#ifndef _AURORA_TWOPOLE_
#define _AURORA_TWOPOLE_
#include "SndBase.h"
#include <cmath>
#include <functional>

namespace Aurora {
enum : int32_t { LP = -1, HP, BP };

template <typename S> inline S id(S s, S dr) {
  (void)dr;
  return s;
}

/** TwoPole class  \n
    2-pole state-variable filter \n
    S: sample type
*/
template <typename S, S (*FN)(S, S) = id> class TwoPole : public SndBase<S> {
  using SndBase<S>::process;
  S Y[2];
  double D[2];
  double W, Fac;
  S ff, dd;
  double piosr;

  S filter(S in, S *y, double *s, double w, double fac, S d, S drv, int32_t typ,
           S m) {
    S lp;
    y[HP] = (in - (d + w) * s[0] - s[1]) * fac;
    S u = w * FN(y[0], drv);
    y[BP] = u + s[0];
    s[0] = y[BP] + u;
    u = w * FN(y[1], drv);
    lp = u + s[1];
    s[1] = lp + u;
    return typ == -1 ? lp * (1 - m) + Y[HP] * m : Y[HP] * (1 - m) + Y[BP] * m;
  }

  void coeffs(S f, S d, S &od, double &w, double &fac, S &of, double ts) {
    w = std::tan(f * ts);
    fac = 1. / (1. + w * d + w * w);
    of = f;
    od = d;
  }

public:
  /** Constructor \n
   sr: sampling rate \n
   vsize: vector size
  */
  TwoPole(S fs = def_sr, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), Y{0}, D{0}, W(0), Fac(0), ff(0), dd(0),
        piosr(M_PI / fs){};

  /** Filter \n
     in: input \n
     f: frequency \n
     d: damping factor (Q reciprocal) \n
     drv: overdrive amount \n
     m: output type (0 - 2: LP(0), HP(1), BP(2))
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S d,
                                   S drv = 0, S m = 0) {
    std::size_t n = 0;
    int32_t typ = m < 1 ? -1 : (m < 2 ? 0 : 1);
    m = m < 0 ? 0 : (m < 1 ? m : (m < 2 ? m - 1 : 1));
    if (f != ff || d != dd)
      coeffs(f, d, dd, W, Fac, ff, piosr);
    auto pf = [&]() {
      return filter(in[n++], Y, D, W, Fac, d, drv + 1, typ, m);
    };
    return process(pf, in.size());
  }

  /** Filter \n
     in: input \n
     f: frequency \n
     d: damping factor (Q reciprocal) \n
     drv: overdrive amount \n
     m: output type (0 - 2: LP(0), HP(1), BP(2))
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f, S d, S drv = 0,
                                   S m = 0) {

    std::size_t n = 0;
    int32_t typ = m < 1 ? -1 : (m < 2 ? 0 : 1);
    m = m < 0 ? 0 : (m < 1 ? m : (m < 2 ? m - 1 : 1));
    auto pf = [&]() {
      if (f[n] != ff || d != dd)
        coeffs(f[n], d, dd, W, Fac, ff, piosr);
      return filter(in[n++], Y, D, W, Fac, d, drv + 1, typ, m);
    };
    return process(pf, in.size() < f.size() ? in.size() : f.size());
  }

  /** reset the filter \n
      fs: sampling rate
   */
  void reset(S fs) {
    piosr = M_PI / fs;
    Y[0] = Y[1] = 0;
    D[0] = D[1] = 0;
    coeffs(ff, dd, dd, W, Fac, ff, piosr);
  }
};
} // namespace Aurora

#endif //_AURORA_TWOPOLE_
