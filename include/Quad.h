// Quad.h:
// Quadrature filter. Outputs an analytic signal.
//
// (c) V Lazzarini, 2022
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

#ifndef _AURORA_QUAD_
#define _AURORA_QUAD_
#include "SndBase.h"
#include <cmath>
#include <complex>

namespace Aurora {

  const double ap1[] = {
    59.017959590337846,
    262.34340692699607,
    1052.8560831644886,
    4223.577583838366,
    17190.389734991037,
    130538.42435798004
  };

  
  const double ap2[] = {
    17.007011830208345,
    129.17600673030512,
    525.775375708461,
    2109.1757722295597,
    8464.591006904155,
    37626.43738022203
  };
  
  /** Quad class  \n
      Quadrature filter \n
      S: sample type
  */
  template <typename S> class Quad : public SndBase<S> {
    using SndBase<S>::get_sig;
    double d1[6], d2[6];
    double c1[6], c2[6];
    std::vector<S> im;
    S ts;

    std::complex<S> filter(S s, double *cr, double *ci, double *dr, double *di) {
      S w;
      std::complex<S> y, x(s, s);
      for(int j = 0; j < 6; j++) {
	w = x.real() + cr[j]*dr[j];
	y.real(dr[j] - cr[j]*w);
	dr[j] = w;
	w = x.imag() + ci[j]*di[j];
	y.imag(di[j] - ci[j]*w);
	di[j] = w;
	x = y;
      }
      return y;  
    }


  public:

    /** Constructor \n
	sr: sampling rate \n
	vsize: vector size
    */   
  Quad(S fs = def_sr, std::size_t vsize = def_vsize)
    : SndBase<S>(vsize), d1{0}, d2{0}, c1{0}, c2{0}, im(vsize), ts(1/fs) 
						       {
							 reset(fs);
						       };


    const std::vector<S> &operator()(const std::vector<S> &in) {
      std::size_t n = 0;
      std::complex<S> cs;
      auto &re = get_sig();
      this->vsize(in.size());
      im.resize(in.size());
      for(auto &s : in) {
	cs = filter(s,c1,c2,d1,d2);
	re[n] = cs.real();
	im[n++] = cs.imag();
      }	
      return re;
    }

    const std::vector<S> &real() const {
      return get_sig();
    }

    const std::vector<S> &imag() const {
      return im;
    }
    
    void reset(S fs) {
      ts = 1/fs;
      S a;
      for(int j = 0; j < 6; j++) {
	a = ap1[j]*ts;
	c1[j] = (1 - a) / (1 + a);
	a = ap2[j]*ts;
	c2[j] = (1 - a) / (1 + a);
	d1[j] = d2[j] = 0;
      }
    }
    
  };
}

#endif
