// SpecShift.h
// Spectral Shifter
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

#ifndef _AURORA_SPECSHIFT_
#define _AURORA_SPECSHIFT_

#include "SpecBase.h"

namespace Aurora {


  /** SpecShift class \n
      Spectral shift \n
      S: sample type
  */    
  template <typename S>
    class SpecShift : public SpecBase<S> {
    using SpecBase<S>::get_spec;
    using SpecBase<S>::fcount_incr;
    Ceps<S> ceps;
    std::vector<S> ftmp;
    S ts;
    bool lock;
      
    const std::vector<specdata<S>> &shift(const std::vector<specdata<S>> &spec,
					  S fscale, S fshift, S forscale,
					  S forshift) {

      auto &buf = get_spec();
      auto size = spec.size();  
      std::size_t n = 0;
      bool preserve;
      S offsr = fshift*size*ts;
   
      std::fill(buf.begin(),buf.end(),Aurora::specdata<S>(0,0));
      if(lock) {
	forscale = 1/fscale;
	forshift = -fshift;
      }
      forshift *= size*ts;
      preserve = forshift != 0 || forscale != 1. ? true : false;
    
      if (preserve) {
	auto &senv = ceps(spec, 30);
	float max = 0;
	for(auto &m : senv) 
	  if(m > max) max = m;
	for(auto &amp : ftmp) {
	  float cf = n/float(size);
	  amp = spec[n].amp();
	  if(senv[n] > 0)	
	    amp *= max/senv[n++];
	}
      }
      	  
      n = 0;	
      for(auto &bin : spec) {
	int k = round(fscale*n + offsr);
	int j = round((1/forscale)*n - forshift);
	auto &senv = ceps.vector();
	if(k > 0  && k < spec.size()) {
	  if(preserve && j > 0  && j < spec.size()
	     && !isnan(senv[j])) {
	      buf[k].amp(ftmp[n]*senv[j]);
	  }
	  else buf[k].amp(bin.amp());
	  buf[k].freq(bin.freq()*fscale + fshift);
	}
	n++;
      }
      return buf;
    }
    

  public:

  SpecShift(S fs = def_sr, std::size_t size = def_fftsize) :
    SpecBase<S>(size), ceps(size), ftmp(size/2 + 1),
      ts(1/fs), lock(false) { };

    const std::vector<specdata<S>> &operator()(const SpecBase<S> &obj,
					       S scl, S shft = 0, S forscl = 0,
					       S forshft = 0) {
      if(obj.framecount() > this->framecount()) {
	fcount_incr();
        return shift(obj.frame(),scl,shft,forscl,forshft);
      } else return get_spec();
    }

    const std::vector<specdata<S>> &operator()(const std::vector<specdata<S>> &spec,
					       S scl, S shft = 0, S forscl = 0,
					       S forshft = 0) {
        return shift(spec,scl,shft,forscl,forshft);
    }


    

    void lock_formants(bool b) { lock = b; }

    void reset(S sr) { ts = 1/sr; }
  };


}
#endif
