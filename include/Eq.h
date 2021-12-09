// Eq.h
// Parametric equaliser
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
#include "SndBase.h"

namespace Aurora {

/** Eq class \n
    Parametric equaliser
    s: sample type
 */
template <typename S> class Eq : SndBase<S> {
  using SndBase<S>::process;
  double z[2];
  double d, a;
  double piosr;
  S ff, bbw;

  S filter(S in, S g, double *z, double d, double a) {
    double w = in + d * (1.0 + a) * z[0] - a * z[1];
    double y = w * a - d * (1.0 + a) * z[0] + z[1];
    z[1] = z[0];
    z[0] = w;
    return (S)(0.5 * (y + in + g * (in - y)));
  }

  void coeffs(S f, S bw) {
    S c = tan(piosr * bw);
    d = std::cos(2 * piosr * f);
    a = (1. - c) / (1. + c);
    ff = f;
    bbw = bw;
  }

public:
  /** Constructor \n
      fs: sampling rate
      vsize: signal vector size
  */
  Eq(S fs = def_sr, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), z{0., 0.}, d(0.), a(0.), piosr(M_PI / fs), ff(0),
        bbw(0){};

  /** Equalisation
     in: input signal
     g: output gain
     fr: centre frequency
     bw: bandwidth
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S g, S fr, S bw) {
    std::size_t n = 0;
    double dd = d, aa = a;
    if (ff != fr || bw != bbw)
      coeffs(fr, bw);
    return process([&]() { return filter(in[n++], g, z, dd, aa); }, in.size());
  }

  /** reset object
      fs: sampling rate
  */
  void reset(S fs) {
    piosr = M_PI / fs;
    coeffs(ff, bbw);
  }
};

} // namespace Aurora
