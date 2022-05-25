// FourPole.h:
// 4-pole lowpass resonating filter
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

#ifndef _AURORA_FOURPOLE_
#define _AURORA_FOURPOLE_

#include "SndBase.h"
#include <cmath>

namespace Aurora {

/** FourPole class  \n
    4-pole lowpass filter \n
    S: sample type
*/
template <typename S> class FourPole : public SndBase<S> {
  using SndBase<S>::process;
  double D[4];
  double A, G[4];
  S ff;
  double piosr;

  S filter(S s, double *d, double *g, double a, S k) {
    S o, u, w;
    S ss = d[3];
    for (int j = 0; j < 3; j++)
      ss += d[j] * g[2 - j];
    o = (g[3] * s + ss) / (1 + k * g[3]);
    u = g[0] * (s - k * o);
    for (int j = 0; j < 3; j++) {
      w = d[j] + u;
      d[j] = u - a * w;
      u = g[0] * w;
    }
    d[3] = g[0] * w - a * o;
    return o;
  }

  void coeffs(S f, double *g, double &a, S &of, double ts) {
    S w = std::tan(f * ts);
    g[0] = w / (1 + w);
    a = (w - 1) / (1 + w);
    g[1] = g[0] * g[0];
    g[2] = g[0] * g[1];
    g[3] = g[0] * g[2];
    of = f;
  }

public:
  /** Constructor \n
   sr: sampling rate \n
   vsize: vector size
  */
  FourPole(S sr = def_sr, int vsize = def_vsize)
      : SndBase<S>(vsize), D{0}, A(0), G{0}, ff(0), piosr(M_PI / sr){};

  /** Filter \n
     in: input \n
     f: cutoff frequency \n
     r: resonance (0-1)
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S r) {
    std::size_t n = 0;
    if (f != ff)
      coeffs(f, G, A, ff, piosr);
    auto pf = [&]() { return filter(in[n++], D, G, A, r * 4); };
    return process(pf, in.size());
  }

  /** Filter \n
   in: input \n
   f: cutoff frequency \n
   r: resonance (0-1)
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &f, S r) {
    std::size_t n = 0;
    auto pf = [&]() {
      if (f[n] != ff)
        coeffs(f[n], G, A, ff, piosr);
      return filter(in[n++], D, G, A, r * 4);
    };
    return process(pf, in.size() < f.size() ? in.size() : f.size());
  }

  /** reset the filter \n
    fs: sampling rate
 */
  void reset(S fs) {
    piosr = M_PI / fs;
    D[0] = D[1] = D[2] = D[3];
    coeffs(ff, G, A, ff, piosr);
  }
};
} // namespace Aurora

#endif // _AURORA_FOURPOLE_
