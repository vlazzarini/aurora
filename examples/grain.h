// grain.h:
// grain processing classes
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

#include "Env.h"
#include "Osc.h"

namespace Aurora {
/** Grain Class: models a grain with a Hanning window envelope */
template <typename S> struct Grain {
  inline static std::vector<S> win{std::vector<S>(0)};
  const std::vector<S> &wave;
  Osc<S, lookupi> env;
  Osc<S, lookupi> osc;
  std::size_t gdr;
  double fs;
  std::size_t t;

  Grain(const std::vector<S> &wve, S sr = def_sr, std::size_t vsize = def_vsize)
      : wave(wve), env(&win, sr, vsize), osc(&wve, sr, vsize), gdr(0), fs(sr),
        t(0) {
    if (win.size() == 0) {
      std::size_t n = 0;
      win.resize(def_ftlen);
      for (auto &s : win)
        s = 0.5 - 0.5 * cos<S>((1. / win.size()) * n++);
    }
  }

  /** trigger a grain of duration d and tab position p (secs) */
  void trigger(S d, S p) {
    gdr = d * fs;
    double ph = p * fs / wave.size();
    osc.phase(ph);
    env.phase(0);
    t = 0;
  }

  /** play grain for set duration with amp a and pitch p */
  auto &operator()(S a, S p) {
    if (t < gdr) {
      S fr = p * fs / wave.size();
      t += osc.vsize();
      return env(osc(a, fr), fs / gdr);
    } else
      return env.clear();
  }
};

/** GrainGen Class: generates streams of grains */
template <typename S> struct GrainGen {
  static S add(S a, S b) { return a + b; }
  std::vector<Grain<S>> slots;
  std::vector<S> mix;
  std::size_t st;
  std::size_t num;
  std::size_t dm;

  GrainGen(const std::vector<S> &wave, std::size_t streams = 16, S sr = def_sr,
           std::size_t decim = def_vsize, std::size_t vsize = def_vsize)
    : slots(streams ? streams : 1, Grain<S>(wave, sr, decim)), mix(vsize), st(0),
      num(0), dm(decim){};

  /** play streams of grains, with amp a, pitch p, grain dur gd (sec),
      density dens (g/sec), and table pos gp (sec) */
  auto &operator()(S a, S p, S dens, S gd, S gp = 0) {
    std::size_t tt = slots[0].fs / dens;
    std::size_t vs = mix.size();
    auto &s = mix;
    for (std::size_t n = 0; n < vs; n+=dm) {
      if (st >= tt) {
        st -= tt;
        slots[num].trigger(gd, gp);
        num = num == slots.size() - 1 ? 0 : num + 1;
      }
      std::fill(s.begin()+n,s.begin()+n+dm,0);
      for (auto &slot : slots) {
        std::size_t j = n;
        for (auto &o : slot(a,p)) 
          s[j++] += o;
       }
      st+=dm;
    }
    return s;
  }
};
} // namespace Aurora
