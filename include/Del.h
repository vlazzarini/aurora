// Del.h
// Generic delay line and associated functions
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

#ifndef _AURORA_DEL_
#define _AURORA_DEL_
#include "SndBase.h"
#include <functional>
#include <iostream>

namespace Aurora {

/** Fixed delay function for Del \n
    S: sample type \n
    nop: no op \n
    wp: reading position (no bounds check) \n
    d: delay line \n
    no1p: no op \n
    returns a sample from the delay line
*/
template <typename S>
inline S fixed_delay(S nop, std::size_t wp, const std::vector<S> &d,
                     std::vector<S> *nop1) {
  (void)nop;
  (void)nop1;
  return d[wp];
}

/* \cond */
template <typename S> inline S rpos(S rp, std::size_t wp, std::size_t ds) {
  rp = wp - rp;
  while (rp < 0)
    rp += ds;
  while (rp >= ds)
    rp -= ds;
  return rp;
}
/* \endcond */

/** Truncating delay function for Del \n
    S: sample type \n
    rp: reading position \n
    wp: write position \n
    d: delay line \n
    nop: no op \n
    returns a sample from the delay line floor(rp) samples behind wp
*/
template <typename S>
inline S vdelay(S rp, std::size_t wp, const std::vector<S> &del,
                std::vector<S> *nop) {
  (void)nop;
  std::size_t ds = del.size();
  return del[(std::size_t)rpos(rp, wp, ds)];
}

/** Interpolation delay function for Del \n
    S: sample type \n
    rp: reading position \n
    wp: write position \n
    d: delay line \n
    nop: no op \n
    returns a sample from the delay line rp samples behind wp, \n
    linearly interpolated
*/
template <typename S>
inline S vdelayi(S rp, std::size_t wp, const std::vector<S> &del,
                 std::vector<S> *nop) {
  (void)nop;
  std::size_t ds = del.size();
  return linear_interp(rpos(rp, wp, ds), del);
}

/** CubicInterpolation delay function for Del \n
    S: sample type \n
    rp: reading position \n
    wp: write position \n
    del: delay line \n
    nop: no op \n
    returns a sample from the delay line rp samples behind wp, \n
    cubic interpolated
*/
template <typename S>
inline S vdelayc(S rp, std::size_t wp, const std::vector<S> &del,
                 std::vector<S> *nop) {
  (void)nop;
  std::size_t ds = del.size();
  return cubic_interp(rpos(rp, wp, ds), del);
}

/** Lowpass-filtered fixed delay function for Del  \n
    nop: no-op \n
    wp: write position \n
    del: delay line \n
    mem: a vector of size 2 with the lp filter state (pos 0) and coeff (pos 1)
   \n returns a convolution sample
*/
template <typename S>
inline S lp_delay(S nop, std::size_t wp, const std::vector<S> &d,
                  std::vector<S> *mem) {
  S ym1 = (*mem)[0];
  S coef = (*mem)[1];
  S x = d[wp];
  S y = (1 + coef) * x - coef * ym1;
  (*mem)[0] = y;
  return y;
}

/** FIR/convolution function for Del  \n
    nop: no-op \n
    wp: write position \n
    del: delay line \n
    ir: impulse response \n
    returns a convolution sample
*/
template <typename S>
inline S fir(S nop, std::size_t wp, const std::vector<S> &del,
             std::vector<S> *ir) {
  auto dl = del.begin() + wp;
  S mx = 0;
  for (auto irs = ir->rbegin(); irs != ir->rend(); irs++, dl++) {
    if (dl == del.end())
      dl = del.begin();
    mx += *dl * *irs;
  }
  return mx;
};

/** Del class \n
    Generic templated delay line \n
    S: sample type
*/
template <typename S, S (*FN)(S, std::size_t, const std::vector<S> &,
                              std::vector<S> *) = fixed_delay>
class Del : public SndBase<S> {
  using SndBase<S>::process;
  S fs;
  std::size_t wp;
  std::vector<S> del;

  S delay(S in, S dt, S fdb, S fwd, std::size_t &p, std::vector<S> *mem) {
    S s = FN(dt, p, del, mem);
    S w = in + s * fdb;
    del[p] = w;
    p = p < del.size() - 1 ? p + 1 : 0;
    return w * fwd + s;
  }

  

public:
  /** Constructor \n
      maxdt: maximum delay time \n
      f: delay lookup function \n
      sr: sampling rate \n
      vsize: vector size
  */
  Del(S maxdt, S sr = def_sr, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), fs(sr), wp(0),
        del(maxdt * fs < 1 ? 1 : maxdt * fs){};

    /** Constructor \n
      maxdt: maximum delay time in Samples \n
      f: delay lookup function \n
      sr: sampling rate \n
      vsize: vector size
  */
   Del(std::size_t maxdt, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), fs(def_sr), wp(0),
        del(maxdt < 1 ? 1 : maxdt){};


  /** Delay \n
      in: audio \n
  */
  const std::vector<S> &operator()(const std::vector<S> &in) {
    std::size_t n = 0, p = wp;
    auto &s = process([&]() { return delay(in[n++], del.size(), 0, 0, p, nullptr); }, in.size());
    wp = p;
    return s;
  }

  
  /** Delay \n
      in: audio \n
      dt: delay time \n
      fdb: feedback gain \n
      fwd: feedforward gain
      mem: optional aux memory
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S dt, S fdb = 0,
                                   S fwd = 0, std::vector<S> *mem = nullptr) {
    std::size_t n = 0, p = wp;
    auto &s = process([&]() { return delay(in[n++], dt*fs, fdb, fwd, p, mem); },
                      in.size());
    wp = p;
    return s;
  }

  /** Delay \n
      in: audio \n
      dt: delay time \n
      fdb: feedback gain \n
      fwd: feedforward gain \n
      mem: optional aux memory
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &dt, S fdb = 0,
                                   S fwd = 0, std::vector<S> *mem = nullptr) {
    std::size_t n = 0, p = wp;

    auto &s = process(
        [&]() {
          auto s = delay(in[n], dt[n] * fs, fdb, fwd, p, mem);
          n++;
          return s;
        },
        in.size() < dt.size() ? in.size() : dt.size());
    wp = p;
    return s;
  }
  
  /** reset the delayline object \n
      maxdt: max delay time \n
      sr: sampling rate
  */
  void reset(S maxdt, S sr) {
    fs = sr;
    wp = 0;
    del.clear();
    del.resize(fs * maxdt);
  }


  std::size_t write_pos() const {
    return wp;
  }

  const std::vector<S> &delayline() const {
    return del;
  }

  
};

/** Tap class \n
    Taps a delay line \n
    S: sample type
*/
template <typename S, S (*FN)(S, std::size_t, const std::vector<S> &,
                              std::vector<S> *) = vdelay>
class Tap : public SndBase<S> {
   using SndBase<S>::process;
   S fs;

  S tap(S dt, const std::vector<S> &del, std::size_t p) {
    return FN(dt, p, del, nullptr);
  }

 public:
  /** Constructor \n
      sr: sampling rate \n
      vsize: vector size
  */
  Tap(S sr = def_sr, std::size_t vsize = def_vsize)
   : SndBase<S>(vsize), fs(sr) { };
   


  /** Tap \n
      del: delay line object \n
      dt: delay time \n
  */ 
   const std::vector<S> &operator()(const Del<S> &del,
                                   const std::vector<S> &dt) {
     std::size_t n = 0;
     int64_t pp = del.write_pos() - del.vsize() - 1;
     auto &d = del.delayline();
    return process(
        [&]() {
	  pp = pp >= -1 ? pp + 1 : pp + 1 + d.size();
          return tap(dt[n++] * fs, d, pp);
        }, dt.vsize());
  }

   /** Tap \n
      del: delay line object \n
      dt: delay time \n
  */ 
   const std::vector<S> &operator()(const Del<S> &del, S dt) {
     int64_t pp = del.write_pos() - del.vsize() - 1;
     auto &d = del.delayline();
     return process([&]() {
	 pp = pp >= -1 ? pp + 1 : pp + 1 + d.size();
	 return tap(dt*fs,d, pp); },
       del.vsize());
   }

  void reset(S sr) {
    fs = sr;
  }

};
 

} // namespace Aurora
#endif // _AURORA_DEL_
