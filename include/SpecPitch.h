// SpecPitch.h
// Pitch Detection
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

#ifndef _AURORA_SPECPITCH_
#define _AURORA_SPECPITCH_

#include "SpecBase.h"

namespace Aurora {

/** SpecPitch class \n
    spectral pitch tracker\n
    S: sample type
*/     
template <typename S> class SpecPitch {
  std::vector<S> peaks;
  std::vector<S> ifacts;
  S cps;

  S frm(S x, S y) {
    return x/y - (int) (x/y);
  }

  S estimate(const std::vector<specdata<S>> &spec, S thresh) {
    std::size_t nbins = spec.size() - 1, np = 0;
    std::size_t n;
    int pp = 0;
    for(n = 1; n < nbins; n++) {
      S a = spec[n].amp();
      if(a > thresh && spec[n-1].amp() < a && spec[n+1].amp() < a) {
	peaks[np++] = spec[n].freq();
	n++;
	if(np == peaks.size()) break;
      }  
    }
    if(!np) return cps;
    bool testa = false , testb = false;
    S end = peaks[0]/20.;
    n = 1;
    for(auto &ifact : ifacts) {
      if(n > end) break;
      S fc, ppk = peaks[0];
      ifact = 0.;
      fc = peaks[0]/n;
      testa = false;
      for(std::size_t k = 0; k < np; k++) {
	S pk = peaks[k];
	S ff = frm(pk,fc);
	int t1 = round(ppk/fc), t2 = round(pk/fc);
	ifact += (ff > 0.5 ? 1. - ff : ff)/pk;
	if(t1 != t2 && t2 - t1 < 3)
	  testa = true;
	ppk = pk;
      }
      if(n == 1) pp = 1;
      else if(ifact < ifacts[pp-1] || (testa && testb)) {
        if(testa) {
	  pp = n;
	  testb = false;
	} else testb = true;
      }
      n++;
    }
    if(!pp) return cps;
    else {
      S fc = peaks[0]/pp, scps = 0;
      for(std::size_t k = 0; k < np; k++)
	scps += peaks[k]/round(peaks[k]/fc);
      cps = scps/np;
    } 
    return cps;
  }


 public:
  
 /** Constructor \n
     npeaks: number of peaks searched for
 */
 SpecPitch(std::size_t npeaks = def_fftsize/4) :
  peaks(npeaks), ifacts(npeaks), cps(0.) { }

  /** Pitch tracking \n
      spec: spectral frame input \n
      thresh: peak finding threshold (0-1)
      returns estimated fundamental frequency
  */
  S operator()(const std::vector<specdata<S>> &spec, S thresh) {
    return estimate(spec, thresh);
  }
  
  /** 
      returns latest estimated fundamental frequency
  */
 S get_cps() { return cps; }
	  
};

}

#endif
