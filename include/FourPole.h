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


#include <cmath>
#include "SndBase.h"

namespace Aurora {


/** FourPole class  \n
    4-pole lowpass filter
*/  
template<typename S>
class FourPole : public SndBase<S> {
  using SndBase<S>::sig;
  double D[4];
  double A, G[4];
  S ff;
  double piosr;

  S filter(S s, double *d, double *g, double a, S k) {
    S o, u, w;
    S ss = d[3];
    for(int j = 0; j < 3; j++) ss += d[j]*g[2-j];
    o = (g[3]*s + ss)/(1 - k*g[3]);
    u = g[0]*(s - k*o);
    for(int j = 0; j < 3; j++) {
      w = d[j] + u;
      d[j] = u - a*w;
      u = g[0]*w;
    }
    d[3] = g[0]*w - a*o;
    return o;
  }

  void coeffs(S f) {
    S g;
    ff = f;
    g = std::tan(f*piosr);
    G[0] = g/(1+g);
    A = A = (g-1)/(1+g);
    G[1] = G[0]*G[0]; 
    G[2] = G[0]*G[1]; 
    G[3] = G[0]*G[2]; 
  }
    
 public:

 /** Constructor
  sr: sampling rate
  vsize: vector size
 */
 FourPole(S sr = def_sr, int vsize = def_vsize) :
  SndBase<S>(vsize), piosr(M_PI/sr) { };

  /** Filter
     in: input
     f: cutoff frequency
     r: resonance (0-1)
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S f, S r) {
    if(f != ff) coeffs(f);
    double *g = G, *d = D, a = A;
    std::size_t n = 0;
    r *= 4;
    for(auto &s:sig) {
      s = filter(in[n++],d,g,a,r);      
    }
    return sig;
  }

    /** Filter
     in: input
     f: cutoff frequency
     r: resonance (0-1)
  */
  const std::vector<S> &operator()(const std::vector<S> &in, const std::vector<S> &f, S r) {
    double *g = G, *d = D, a = A;
    std::size_t n = 0;
    r *= 4;
    for(auto &s:sig) {
      if(f[n] != ff) coeffs(f[n]);
      s = filter(in[n++],d,g,a,r);      
    }
    return sig;
  }
};
}
