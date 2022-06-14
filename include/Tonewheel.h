// Tonewheek.h:
// Tonewheel model
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

#ifndef _AURORA_TW_
#define _AURORA_TW_

#include "SndBase.h"
#include "Osc.h"

namespace Aurora {
  template <typename S>
    struct Tonewheel : public SndBase<S> {
    static constexpr int64_t maxlen = 0x100000000; // max tab len 2^32
    using SndBase<S>::process;
    const std::vector<S> *tab;
    int64_t fac;
    int lobits; // how many bits are used in indexing
    int lomask; // mask for frac extraction
    S lofac; // fac for frac extract
    
    S lookup(unsigned int phs) {
      auto &t = *tab;
      float frac = (phs & lomask)*lofac; // frac part of index
      unsigned int ndx = phs >> lobits;  // index
      S s = t[ndx] + frac*(t[ndx+1] - t[ndx]);  // lookup
      return s;
    }

  public:
  Tonewheel(const std::vector<S> *w = NULL, S ratio = 1,
	    std::size_t vsize = def_vsize) : SndBase<S>(vsize),
      tab(w), fac(ratio*maxlen), lobits(0) {
      if(w) {
      for(int64_t t = tab->size()-1; (t & maxlen) == 0; t <<= 1) 
	lobits += 1;
      lomask = (1 << lobits) - 1;
      lofac = 1.f/(lomask + 1);
      }
    }

    void set_ratio(S r) { fac = r*maxlen; }
    void set_table(const std::vector<S> *w) {
      tab = w;
      for(int64_t t = tab->size()-1; (t & maxlen) == 0; t <<= 1) 
	lobits += 1;
      lomask = (1 << lobits) - 1;
      lofac = 1.f/(lomask + 1);
    }

    const std::vector<S> &operator() (const std::vector<S> &phs) {
      std::size_t n  = 0;
      return process(
        [&]() {
          return lookup((unsigned int)(phs[n++]*fac));
        },
        phs.size());
    }
  };

template <typename S>
struct SineTab {
  std::vector<S> tab;
  SineTab() : tab(pow(2,14)+1) {
    std::size_t n = 0;
    for(auto &s : tab)
      s = std::sin((n++)*twopi/(tab.size()-1));
  }
};


template <typename S>
struct SqrTab {
  std::vector<S> tab;
  SqrTab() : tab(pow(2,14)+1) {
    std::size_t n = 0;
    for(auto &s : tab) {
      s = std::sin((n++)*twopi/(tab.size()-1));
      s += std::sin(3*(n++)*twopi/(tab.size()-1))/3;
      s += std::sin(5*(n++)*twopi/(tab.size()-1))/5;
    }
    S max = 0;
    for(auto &s : tab) if(fabs(s) > max) max  = fabs(s);
    for(auto &s : tab) s *= 1/max;
  }
}; 

template <typename S>
class Tonegen {
  const static SineTab<S> stab;
  const static SqrTab<S>  sqtab;
  std::vector<Tonewheel<S>> wheels;
  std::vector<Osc<S,phase>> phs;
  S ffs[12];

  void run(float scal, std::size_t vsiz) {
    std::size_t n = 0;
    for(auto &p : phs) {
      p.vsize(vsiz);
      p(scal,ffs[n++]);
    }
    n = 0;
    for(auto &w : wheels) {
      w(phs[n%12].vector());
      n++;
     }
  }

 public:
 Tonegen() : wheels(91), phs(12) {
    for(std::size_t n = 0; n < wheels.size(); n++) {
      if(n < 12) {
	wheels[n].set_table(&(sqtab.tab));
        ffs[n] = pow(2,n/12.)*32.69;
      }
      else wheels[n].set_table(&(stab.tab));
      wheels[n].set_ratio(pow(2,n/12));
    }
 
  }

  void operator()(float scal, std::size_t vsize) { run(vsize); }
  const std::vector<S> &wheel(float scal, std::size_t num) {
    return wheels[num].vector();
  }
 
  };

template <typename S>
const SineTab<S> Tonegen<S>::stab;

template <typename S>
const SqrTab<S> Tonegen<S>::sqtab;
}

#endif
