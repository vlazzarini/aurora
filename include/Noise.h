// Noise.h:
// Noise Generators
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

#ifndef _AURORA_NOISE_
#define _AURORA_NOISE_
#include "SndBase.h"

namespace Aurora {

  template <typename S> S white(S a) {
    return a*(2.*std::rand()/((S)RAND_MAX) - 1.);
  }

/** Noise class  \n
    Generic noise generator \n
    S: sample type \n
    FN: noise function
*/
template <typename S, S (*FN)(S) = white>
class Noise : public SndBase<S> {
  using SndBase<S>::process;
  S fs;
  S incr;
  S ov;
  std::size_t t;

  S sample(S a, std::size_t pp, bool interp) {
    if(++t >= pp)
      {
      S nv = FN(a);
      incr = (nv - ov)/pp;
      t = 0;
      if(!interp) ov = nv;
     }
    ov = interp ? ov + incr : ov;
    return ov;
  
  }


 public:
 Noise(S sr = def_sr,std::size_t vsize = def_vsize) 
  : SndBase<S>(vsize), fs(sr), incr(0), ov(0), t(0) { }

  const std::vector<S> & operator()(S a, S f, bool interp = false) {
    std::size_t pp = fs/(f > 0 ? f : 0.000001);
    return process([&]() {
	return sample(a, pp, interp);
	}, 0);
  }

  const std::vector<S> & operator()(S a) {
    std::size_t tt = t;
    return process([&]() {
	return sample(a,0,0); }, 0);
  }

};
 
}

#endif
