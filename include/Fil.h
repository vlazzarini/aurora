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
template <typename S> inline S reson(S in, double *c, double *d) {
  S y = in * c[0] - d[0] * c[1] - d[1] * c[2];
  d[1] = d[0];
  d[0] = y;
  return y;
}

/** resonator coefficients function for Fil
    no scaling
 */
template <typename S> inline void reson_cfs(S f, S bw, S fs, double *c) {
  c[2] = std::exp(-bw * twopi / fs);
  c[1] = (-4 * c[2] / (1. + c[2])) * std::cos(twopi * f / fs);
}

/** resonator coefficients function for Fil
    scaling type 1
 */
template <typename S> inline void reson_cfs1(S f, S bw, S fs, double *c) {
  reson_cfs(f, bw, fs, c);
  c[0] = (1 - c[2]) * sqrt(1.0 - c[1] * c[1] / (4 * c[2]));
}

/** resonator coefficients function for Fil
    scaling type 2
 */
template <typename S> inline void reson_cfs2(S f, S bw, S fs, double *c) {
  reson_cfs(f, bw, fs, c);
  double rsqp1 = c[2] + 1;
  c[0] = sqrt((rsqp1 * rsqp1 - c[1] * c[1]) * (1 - c[2]) / rsqp1);
}

/** DF-I second-order section
 */
template <typename S> inline S dfI(S in, double *c, double *d) {
  S y = in * c[0] + d[0] * c[1] + d[1] * c[2] - d[2] * c[3] - d[3] * c[4];
  d[1] = d[0];
  d[0] = in;
  d[3] = d[2];
  d[2] = y;
  return y;
}

/** DF-II second-order section
 */
template <typename S> inline S dfII(S in, double *c, double *d) {
  S w = in - d[0] * c[3] - d[1] * c[4];
  S y = w * c[0] + d[0] * c[1] + d[1] * c[2];
  d[1] = d[0];
  d[0] = w;
  return y;
}

/** Second-order lowpass filter coeffs
 */
template <typename S> inline void lp_cfs(S f, S bw, S fs, double *c) {
  double w = 1 / tan(M_PI * f / fs);
  double sqw = sqrt(2.) * w;
  double wsq = w * w;
  c[0] = 1 / (1 + sqw + wsq);
  c[1] = 2 * c[0];
  c[2] = c[0];
  c[3] = 2 * (1. - sqw) * c[0];
  c[4] = (1 - sqw + wsq) * c[0];
}

/** Second-order hipass filter coeffs
 */
template <typename S> inline void hp_cfs(S f, S bw, S fs, double *c) {
  double w = tan(M_PI * f / fs);
  double sqw = sqrt(2.) * w;
  double wsq = w * w;
  c[0] = 1 / (1 + sqw + wsq);
  c[1] = -2 * c[0];
  c[2] = c[0];
  c[3] = 2 * (wsq - 1) * c[0];
  c[4] = (1 - sqw + wsq) * c[0];
}

/** Second-order bandpass filter coeffs
 */
template <typename S> inline void bp_cfs(S f, S bw, S fs, double *c) {
  double w = 1. / tan(M_PI * bw / fs);
  double cosw = 2. * cos(2 * M_PI * f / fs);
  c[0] = 1. / (1 + w);
  c[1] = 0;
  c[2] = -c[0];
  c[3] = -w * cosw * c[0];
  c[4] = (w - 1) * c[0];
}

/** Second-order notch filter coeffs
 */
template <typename S> inline void br_cfs(S f, S bw, S fs, double *c) {
  double w = tan(M_PI * bw / fs);
  double cosw = 2. * cos(2 * M_PI * f / fs);
  c[0] = 1 / (1 + w);
  c[1] = -cosw * c[0];
  c[2] = c[0];
  c[3] = c[1];
  c[4] = (1 - w) * c[0];
}

/** Fil class  \n
    Generic filter (first or second-order) \n
    S: sample type
*/
template <typename S, void (*CF)(S, S, S, double *),
          S (*FN)(S, double *, double *)>
class Fil : public SndBase<S> {
  using SndBase<S>::process;
  double d[4];
  double c[5];
  S ff;
  S bbw;
  S fs;

  S filter(S s, double *c, double *d) { return FN(s, c, d); }

  void coeffs(S f, S bw, S fs, double *c) {
    CF(f, bw, fs, c);
    ff = f;
    bbw = bw;
  }

public:
  /** Constructor \n
   fs: sampling rate \n
   vsize: vector size
  */
  Fil(S ffs = def_sr, int vsize = def_vsize)
      : SndBase<S>(vsize), d{0}, c{0}, ff(0), bbw(0), fs(ffs){};

  /** Filter \n
     in: input \n
     f: cutoff frequency \n
     bw: bandwidth \n;
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S bw = 0) {
    std::size_t n = 0;
    double *D = d, *C = c;
    if (f != ff || bw != bbw)
      coeffs(f, bw, fs, c);
    auto &s = process([&]() { return filter(in[n++], C, D); }, in.size());
    return s;
  }

  /** Filter \n
   in: input \n
   f: cutoff frequency \n
   bw: bandwidth \n;
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f, S bw = 0) {
    std::size_t n = 0;
    double *D = d, *C = c;
    auto &s = process(
        [&]() {
          if (f[n] != ff)
            coeffs(f, bw, fs, C);
          return filter(in[n++], C, D);
        },
        in.size() < f.size() ? in.size() : f.size());
    return s;
  }

  /** reset the filter \n
       fs: sampling rate
    */
  void reset(S ffs) {
    d[0] = d[1] = d[2] = d[3] = 0;
    coeffs(ff, bbw, fs, c, d);
    fs = ffs;
  }
};
} // namespace Aurora

#endif // _AURORA_ONEPOLE_
