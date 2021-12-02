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
  rp = wp - rp;
  while (rp < 0)
    rp += ds;
  while (rp >= ds)
    rp -= ds;
  return del[(std::size_t)rp];
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
  rp = wp - rp;
  while (rp < 0)
    rp += ds;
  while (rp >= ds)
    rp -= ds;
  return linear_interp(rp, del);
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
  std::function<S(S, std::size_t &, const std::vector<S> &)> fun;

  S delay(S in, S dt, S fdb, S fwd, std::size_t &p) {
    S s = fun(dt * fs, p, del);
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
  Del(S maxdt, std::function<S(S, std::size_t, const std::vector<S> &)> f,
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
    const std::vector<S> &s =
        process([&]() { return delay(in[n++], dt, fdb, fwd, p); }, in.size());
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
    const std::vector<S> &s = process(
        [&]() {
          auto s = delay(in[n], dt[n], fdb, fwd, p);
          n++;
          return s;
        },
        in.size() < dt.size() ? in.size() : dt.size());
    wp = p;
    return s;
  }
};

} // namespace Aurora
#endif // _AURORA_DEL_
