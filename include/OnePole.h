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
  using SndBase<S>::process;
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

  void coeffs(S f, double &g, double &a, S &of, double ts) {
    S w = std::tan(f * ts);
    g = w / (1 + w);
    a = (w - 1) / (1 + w);
    of = f;
  }

public:
  /** Constructor \n
   fs: sampling rate \n
   vsize: vector size
  */
  OnePole(S fs = def_sr, int vsize = def_vsize)
      : SndBase<S>(vsize), D(0), A(0), G(0), ff(0), piosr(M_PI / fs){};

  /** Filter \n
     in: input \n
     f: cutoff frequency \n
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f) {
    double d = D;
    std::size_t n = 0;
    if (f != ff)
      coeffs(f, G, A, ff, piosr);
    auto &s = process([&]() { return filter(in[n++], d, G, A); }, in.size());
    D = d;
    return s;
  }

  /** Filter \n
   in: input \n
   f: cutoff frequency \n
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f) {
    double d = D;
    std::size_t n = 0;
    auto &s = process(
        [&]() {
          if (f[n] != ff)
            coeffs(f[n], G, A, ff, piosr);
          return filter(in[n++], d, G, A);
        },
        in.size() < f.size() ? in.size() : f.size());
    D = d;
    return s;
  }

  /** reset the filter \n
       fs: sampling rate
    */
  void reset(S fs) {
    D = 0;
    piosr = M_PI / fs;
    coeffs(ff, G, A, ff, piosr);
  }
};
} // namespace Aurora

#endif // _AURORA_ONEPOLE_
