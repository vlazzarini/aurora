// OnePole.h:
// First-order lowpass filter
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

/** OnePole class  \n
    First-order lowpass filter \n
    S: sample type
*/
template <typename S> class OnePole : public SndBase<S> {
  using SndBase<S>::sig;
  double D;
  double A, G;
  S ff;
  double piosr;

  S filter(S s, double &d, double g, double a) {
    S u = g * s;
    S y = u + d;
    d = u - a * y;
    return y;
  }

  void coeffs(S f) {
    S g = std::tan(f * piosr);
    G = g / (1 + g);
    A = (g - 1) / (1 + g);
    ff = f;
  }

public:
  /** Constructor \n
   sr: sampling rate \n
   vsize: vector size
  */
  OnePole(S sr = def_sr, int vsize = def_vsize)
      : SndBase<S>(vsize), D(0), A(0), G(0), ff(0), piosr(M_PI / sr){};

  /** Filter \n
     in: input \n
     f: cutoff frequency \n
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f) {
    if (f != ff)
      coeffs(f);
    double d = D, a = A, g = G;
    std::size_t n = 0;
    for (auto &s : sig) {
      s = filter(in[n++], d, g, a);
    }
    D = d;
    return sig;
  }

  /** Filter \n
   in: input \n
   f: cutoff frequency \n
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f) {
    double d = D, &a = A, &g = G;
    std::size_t n = 0;
    for (auto &s : sig) {
      if (f[n] != ff)
        coeffs(f[n]);
      s = filter(in[n++], d, g, a);
    }
    D = d;
    return sig;
  }
};
} // namespace Aurora

#endif // _AURORA_ONEPOLE_