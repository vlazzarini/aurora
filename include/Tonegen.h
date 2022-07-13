// Tonegen.h:
// Tone generator models
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
#include "BlOsc.h"

namespace Aurora {
  template <typename S>
    struct Lookup : public SndBase<S> {
    static constexpr int32_t maxlen = 0x40000000; // max tab len 2^30
    static constexpr int32_t phmsk = maxlen - 1;
    using SndBase<S>::process;
    const std::vector<S> *tab;
    int64_t fac;
    int lobits; // how many bits are used in indexing
    int lomask; // mask for frac extraction
    S lofac; // fac for frac extract
    
    S lookup(int64_t phs) {
      auto &t = *tab;
      float frac = (phs & lomask)*lofac; // frac part of index
      int64_t ndx = (phs & phmsk) >> lobits;  // index
      S s = t[ndx] + frac*(t[ndx+1] - t[ndx]);  // lookup
      return s;
    }

  public:
  Lookup(const std::vector<S> *w = NULL, S ratio = 1,
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
      lobits = 0;
      for(int64_t t = tab->size()-1; (t & maxlen) == 0; t <<= 1) 
	lobits += 1;
      lomask = (1 << lobits) - 1;
      lofac = 1.f/(lomask + 1);
    }


    void swap_table(const std::vector<S> *w) { tab = w; }
      
    const std::vector<S> &operator() (const std::vector<S> &phs) {
      std::size_t n  = 0;
      return process(
        [&]() {
          return lookup((phs[n++]*fac));
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
      s = std::sin((n)*twopi/(tab.size()-1));
      s += std::sin(3*(n)*twopi/(tab.size()-1))/3;
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
  std::vector<Lookup<S>> wheels;
  std::vector<Osc<S,phase>> phs;
  S ffs[12];

  void run(std::size_t vsiz) {
    std::size_t n = 0;
    for(auto &p : phs) {
      p.vsize(vsiz);
      p(1,ffs[n++]);
    }
    n = 0;
    for(auto &w : wheels) {
      w(phs[n%12].vector());
      n++;
     }
  }

 public:
 Tonegen() : wheels(91), phs(12), ffs{0.817307692,0.865853659,0.917808219,
       0.972222222,1.03,1.090909091,1.15625,1.225,1.297297297,1.375,
       1.456521739,1.542857143}
       {
    for(std::size_t n = 0; n < wheels.size(); n++) {
      if(n < 12) {
	wheels[n].set_table(&(sqtab.tab));
        ffs[n] *= 40; // 2 teeth @ 20 rev/s 
      }
      else wheels[n].set_table(&(stab.tab));
      wheels[n].set_ratio(pow(2,n/12));
    }
 
  }

  void operator()(std::size_t vsize) { run(vsize); }
  const std::vector<S> &wheel(std::size_t num) {
    return wheels[num].vector();
  }

  void reset(S fs) {
    for(auto &p : phs) p.reset(fs);
  }
 
  };

template <typename S>
const SineTab<S> Tonegen<S>::stab;

template <typename S>
const SqrTab<S> Tonegen<S>::sqtab;
 
  template <typename S>
   class Wavegen  {
    TableSet<S> waveset;
    std::vector<Osc<S,phase>> phs;
    Lookup<S> tread;
    std::vector<double> freq;
    std::vector<S> mix;

    void run(std::size_t vsiz, S detun = 1.f) {
      std::size_t n = 0;
      for(auto &p : phs) {
	p.vsize(vsiz);
	p(1,freq[n++]*detun);
      }
    }

  public:

  Wavegen(int type = SAW) : waveset(type), phs(12), freq(128), mix(def_vsize) {
       std::size_t n = 0;
       float base = 8.175798915643707;
       for(auto &f : freq) 
	 f = base*pow(2,n++/12.);
        waveset.guardpoint();
	tread.set_table(&waveset.func(freq[0]));
    }

    const std::vector<S> &operator()(std::size_t vsize, S detun = 1.f) {
      mix.resize(vsize);
      tread.vsize(vsize);
      std::fill(mix.begin(), mix.end(), 0);
      run(vsize,detun);
      return mix;
    }
				   
    const std::vector<S> &tone(std::size_t note, S amp) {
      if(note < 128) {
      tread.set_ratio(freq[note]/freq[note%12]);
      tread.swap_table(&waveset.func(freq[note]));      
      auto &tmp = tread(phs[note%12].vector());
      std::size_t n = 0;
      for(auto &s : mix) s += tmp[n++]*amp;
      return tmp;
      } else return mix;
    }

    void reset(S fs, int type = SAW) {
      waveset.reset(type, fs);
           waveset.guardpoint();   
      tread.set_table(&waveset.func(freq[0]));
      for(auto &p : phs) p.reset(fs);
    } 
  };


}

#endif
