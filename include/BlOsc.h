// BlOsc.h
// Wavetable sets
// Bandlimited Wavetable Oscillator
//
// (c) V Lazzarini, 2021
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors
// may be used to endorse or promote products derived from this software without
// specific prior written permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

#ifndef _AURORA_BLOSC_
#define _AURORA_BLOSC_

#include "Osc.h"

namespace Aurora {
  enum {SAW = 0, SQUARE, TRIANGLE, PULSE};
  const int octs = 10;
  const double base = 20;
  const int def_ftlen = 16384;

  /** TableSet class \n
      Creates a set of tables for BlOsc.
  */
  template<typename S> class TableSet {
    std::size_t tlen;
    std::vector<std::vector<S>> waves;

    void norm(std::vector<S> &wave) {
      S max = 0.;
      for (auto s: wave)
        if (s > max) max = s;
      for (auto &s: wave)
        s *= 1./max;
    }
  
    void fourier(std::vector<S> &wave, std::size_t nh, uint32_t type) {
      double ol = 1./tlen;
      std::vector<S> h(nh); 
      std::iota(h.begin(), h.end(), 1);
      std::function<S(S)> fun = type < TRIANGLE ? sin<S> : cos<S>;
      for (auto &n: h) {
        std::size_t k = 0;
        for (auto &s: wave) {
          S a;
          switch(type) {
          case SAW:
            a = 1./n;
            break;
          case SQUARE:
            a = (int) n % 2 ? 1./n : 0.;
            break;
          case TRIANGLE:
            a = (int) n % 2 ? 1./(n*n) : 0.;
            break;
          default:
            a = 1.;
          }
          s += a*fun(n*ol*k++);
        }
      }
      norm(wave);
    }

  public:
    /** Constructor \n
        type: wave type (SAW, SQUARE, TRIANGLE, PULSE) \n
        fs: sampling rate for which these will be built \n
        len: table length
    */  
  TableSet(uint32_t type, S fs = def_sr, std::size_t len = def_ftlen)
    : tlen(len), waves(octs,std::vector<S>(len)) {
      S nyq = fs/2.;
      int i = 0;
      for (auto &wave:waves) {
        double fr = base*pow(2, i++);
        fourier(wave,nyq/fr,type);
      }
    }

    /** Table selection \n
        f: fundamental frequency used for playback
    */
    const std::vector<S> &select(S f) const {
      int32_t num = (int32_t) round(std::log2(f/base));
      return num < waves.size() ? waves[num] : waves.back();
    } 
  };

  /** BlOsc class \n
      Bandlimited wavetable oscillator.
  */
  template<typename S> class BlOsc : public Osc<S> {
    using Osc<S>::ph;
    using Osc<S>::ts;
    const TableSet<S> &tset;
    
    virtual S synth(S a,S f,double &phs) {
      const std::vector<S> &t = tset.select(f);
      size_t len = t.size();
      size_t posi = (size_t) phs;
      double frac = phs - posi;
      S s = a*(t[posi] + frac*((posi != len ? t[posi+1] : t[posi]) - t[posi]));
      phs += f*ts*len;
      while (phs < 0) phs += len; 
      while (phs >= len) phs -= len;
      return s;
    }
       
  public:
    /** Constructor \n
        t: wavetable set \n
        fs: sampling rate \n
        vsize: vector size
    */
  BlOsc(const TableSet<S>& t, S fs = (S) def_sr,
        std::size_t vsize = def_vsize) :
    Osc<S>(nullptr,fs,vsize), tset(t) { };
  };
}

#endif // _AURORA_BLOSC_
