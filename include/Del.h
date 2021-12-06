// Del.h
// Generic delay line and associated functions
//
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
    rp: no op \n
    wp: reading position (no bounds check) \n
    d: delay line \n
    returns a sample from the delay line
*/
template <typename S>
S fixed_delay(S rp, std::size_t wp, const std::vector<S> &d) {
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
    returns a sample from the delay line floor(rp) samples behind wp
*/
template <typename S>
S vdelay(S rp, std::size_t wp, const std::vector<S> &del) {
  std::size_t ds = del.size();
  return del[(std::size_t)rpos(rp, wp, ds)];
}

/** Interpolation delay function for Del \n
    S: sample type \n
    rp: reading position \n
    wp: write position \n
    d: delay line \n
    returns a sample from the delay line rp samples behind wp, \n
    linearly interpolated
*/
template <typename S>
S vdelayi(S rp, std::size_t wp, const std::vector<S> &del) {
  std::size_t ds = del.size();
  return linear_interp(rpos(rp, wp, ds), del);
}

/** CubicInterpolation delay function for Del \n
    S: sample type \n
    rp: reading position \n
    wp: write position \n
    d: delay line \n
    returns a sample from the delay line rp samples behind wp, \n
    cubic interpolated
*/
template <typename S>
S vdelayc(S rp, std::size_t wp, const std::vector<S> &del) {
  std::size_t ds = del.size();
  return cubic_interp(rpos(rp, wp, ds), del);
}

/** Generating function for lpf delay function  \n
    s: externally-defined filter state \n
    c: lp filter coef [c = sqrt(a*a - 1) - a, with a  = 2 - cos(w)] \n
    f: delay function
    returns a delay function for use in Del
*/
template <typename S>
std::function<S(S, std::size_t, const std::vector<S> &)>
lpdelay_gen(S &s, double &c,
            std::function<S(S, std::size_t, const std::vector<S> &)> f) {
  return [&](S rp, std::size_t wp, const std::vector<S> &d) -> S {
    auto o = f(rp, wp, d);
    return (s = o * (1 + c) - s * c);
  };
}

/** Generating function for FIR/convolution function  \n
    ir: impulse response \n
    returns a delay function for use in Del
*/
template <typename S>
std::function<S(S, std::size_t, const std::vector<S> &)>
fir_gen(const std::vector<S> &ir) {
  return [&](S rp, std::size_t wp, const std::vector<S> &del) -> S {
    auto dl = del.begin() + wp;
    S mx = 0;
    for (auto irs = ir.rbegin(); irs != ir.rend(); irs++, dl++) {
      if (dl == del.end())
        dl = del.begin();
      mx += *dl * *irs;
    }
    return mx;
  };
}

/** Del class \n
    Generic delay line \n
    S: sample type
*/
template <typename S> class Del : SndBase<S> {
  using SndBase<S>::process;
  S fs;
  std::size_t wp;
  std::vector<S> del;
  std::function<S(S, std::size_t, const std::vector<S> &)> fun;

  S delay(S in, S dt, S fdb, S fwd, std::size_t &p,
          std::function<S(S, std::size_t, const std::vector<S> &)> fn) {
    S s = fn(dt * fs, p, del);
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
  Del(S maxdt, const std::function<S(S, std::size_t, const std::vector<S> &)> f,
      S sr = def_sr, std::size_t vsize = def_vsize)
      : SndBase<S>(vsize), fs(sr), wp(0), del(maxdt * fs), fun(f){};

  /** Constructor \n
      dt: fixed delay time \n
      f: delay lookup function \n
      sr: sampling rate \n
      vsize: vector size
  */
  Del(S dt, S sr = def_sr, std::size_t vsize = def_vsize)
      : Del(dt, fixed_delay<S>, sr, vsize){};

  /** Delay \n
      in: audio \n
      dt: delay time \n
      fdb: feedback gain \n
      fwd: feedforward gain
  */
  const std::vector<S> &operator()(const std::vector<S> &in, S dt, S fdb = 0,
                                   S fwd = 0) {
    std::size_t n = 0, p = wp;
    auto &s = process([&]() { return delay(in[n++], dt, fdb, fwd, p, fun); },
                      in.size());
    wp = p;
    return s;
  }

  /** Delay \n
      in: audio \n
      dt: delay time \n
      fdb: feedback gain \n
      fwd: feedforward gain
  */
  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &dt, S fdb = 0,
                                   S fwd = 0) {
    std::size_t n = 0, p = wp;
    auto &s = process(
        [&]() {
          auto s = delay(in[n], dt[n], fdb, fwd, p, fun);
          n++;
          return s;
        },
        in.size() < dt.size() ? in.size() : dt.size());
    wp = p;
    return s;
  }

  /** set the delay function \n
      f: delay function to be used
  */
  void func(const std::function<S(S, std::size_t, const std::vector<S> &)> f) {
    fun = f;
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
};

} // namespace Aurora
#endif // _AURORA_DEL_
