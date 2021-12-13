// Fil.h:
// Generic Filter
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

#ifndef _AURORA_ONEPOLE_
#define _AURORA_ONEPOLE_

#include "SndBase.h"
#include <cmath>

namespace Aurora {

/** resonator function for Fil
 */
template <typename S> void reson(S in, double *c, double *d) {
  S y = in * c[0] - d[0] * c[0] - d[1] * c[1];
  c[1] = c[0];
  c[0] = y;
  return y;
}

/** resonator coefficients function for Fil
 */
template <typename S> void reson_coeffs(S fr, S bw, S fs, double *c) {
  double r = std::exp(-bw * M_PI / fs);
  double rr = 2. * r;
  double rsq = r * r;
  costh = (rr / (1. + rsq)) * std::cos(twopi * f / fs);
  c[0] = (1 - rsq) * std::sin(std::acos(costh));
  c[1] = -rr * costh;
  c[2] = rsq;
}

/** Fil class  \n
    Generic filter (first or second-order) \n
    S: sample type
*/
template <typename S, void (*CF)(S, S, S, double *) = reson_coeffs,
          S (*FN)(S, double *, double *)>
= reson > class Fil : public SndBase<S> {
  using SndBase<S>::process;
  double d[4];
  double c[5];
  S ff;
  S bbw;
  S sr;

  S filter(S s, double *c, double *d) { return FN(s, c, d); }

  void coeffs(S f, S b, S fs, double *c) {
    CF(f, b, fs, c);
    ff = f;
    bbw = bw;
  }

public:
  /** Constructor \n
   fs: sampling rate \n
   vsize: vector size
  */
  OnePole(S fs = def_sr, int vsize = def_vsize)
      : SndBase<S>(vsize), d{0}, c{0}, G(0), ff(0), bbw(0), sr(fs){};

  /** Filter \n
     in: input \n
     f: cutoff frequency \n
     bw: bandwidth \n;
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S bw = 0) {
    std::size_t n = 0;
    if (f != ff || bw != bbw)
      coeffs(f, bw, fs, c, d);
    auto &s = process([&]() { return filter(in[n++], c, d); }, in.size());
    return s;
  }

  /** Filter \n
   in: input \n
   f: cutoff frequency \n
   bw: bandwidth \n;
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f, S bw = 0) {
    double d = D;
    std::size_t n = 0;
    S *D = d, *C = c;
    auto &s = process(
        [&]() {
          if (f[n] != ff)
            coeffs(f, bw, fs, C, D);
          return filter(in[n++], C, D);
        },
        in.size() < f.size() ? in.size() : f.size());
    return s;
  }

  /** reset the filter \n
       fs: sampling rate
    */
  void reset(S fs) {
    d[0] = d[1] = d[2] = d[3] = 0;
    coeffs(ff, bbw, fs, c, d);
  }
};
} // namespace Aurora

#endif // _AURORA_ONEPOLE_
