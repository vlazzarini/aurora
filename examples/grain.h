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
  bool off;

  Grain(const std::vector<S> &wve, S sr = def_sr, std::size_t vsize = def_vsize)
      : wave(wve), env(&win, sr, vsize), osc(&wve, sr, vsize), gdr(0), fs(sr),
       t(0), off(false) {
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
    off = false;
  }

  void reset(S sr) {
    fs = sr;
    osc.reset(fs);
    env.reset(fs);
  }

  void vsize(std::size_t vs) {
    osc.vsize(vs);
    env.vsize(vs);
  }

  /** play grain for set duration with amp a and pitch p */
  auto &operator()(S a, S p) {
    if (t < gdr) {
      S fr = p * fs / wave.size();
      t += osc.vsize();
      return env(osc(a, fr), fs / gdr);
    } else {
      if(!off) {
        env.clear();
        off = true;
      }
      return env.vector();
    }
  }

  /** play grain for set duration with AM and FM */
  auto &operator()(const std::vector<S> am, const std::vector<S> fm) {
    if (t < gdr) {
      t += osc.vsize();
      return env(osc(am, fm), fs / gdr);
    } else {
      if(!off) {
        env.clear();
        off = true;
      }
      return env.vector();
    }
  }

    /** play grain for set duration with AM and PM */
  auto &operator()(const std::vector<S> a, S f, const std::vector<S> pm) {
    if (t < gdr) {
      t += osc.vsize();
      return env(osc(a, f, pm), fs / gdr);
    } else {
      if(!off) {
        env.clear();
        off = true;
      }
      return env.vector();
    }
  }

  auto &operator()(S a, S f, S pm) {
    if (t < gdr) {
      t += osc.vsize();
      return env(osc(a, f), fs / gdr);
    } else {
      if(!off) {
        env.clear();
        off = true;
      }
      return env.vector();
    }
  }

  
};

/** GrainGen Class: generates streams of grains */
template <typename S> struct GrainGen {
  static S add(S a, S b) { return a + b; }
  std::vector<Grain<S>> slots;
  std::vector<S> mixl;
  std::vector<S> mixr;
  std::size_t st;
  std::size_t num;
  std::size_t dm;
  std::size_t dmr;

  GrainGen(const std::vector<S> &wave, std::size_t streams = 16, S sr = def_sr,
           std::size_t decim = def_vsize, std::size_t vsize = def_vsize)
  : slots(streams ? streams : 1, Grain<S>(wave, sr, decim)), mixl(vsize), mixr(vsize), st(0),
    num(0), dm(decim), dmr(dm/vsize) {};

  void reset(S fs){
    for (auto &grain: slots) grain.reset(fs);
  }

  /** play streams of grains, with amp a, pitch p, grain dur gd (sec),
      density dens (g/sec), and table pos gp (sec) */
  auto &operator()(S a, S p, S dens, S gd, S gp = 0, std::size_t vs = def_vsize) {
    auto &grains = slots;
    auto &s = mixl;
    std::size_t ddm = dmr*vs;
    std::size_t tt = grains[0].fs / dens;
    s.resize(vs);
    for (std::size_t n = 0; n < vs; n+=ddm) {
      if (st >= tt) {
        st -= tt;
        grains[num].trigger(gd, gp);
        num = num == slots.size() - 1 ? 0 : num + 1;
      }
      std::fill(s.begin()+n,s.begin()+n+ddm,0);
      for (auto &grain: grains) {
        std::size_t j = n;
	grain.vsize(ddm);
        for (auto &o : grain(a,p)) 
          s[j++] += o;
       }
      st+=ddm;
    }
    return s;
  }


  /** play streams of grains, with amp am, freq f, pm pm, grain dur gd (sec),
      density dens (g/sec), and table pos gp (sec) */
  auto &operator()(const std::vector<S> am, S f, const std::vector<S> pm, S pan, S dens, S gd,
		   S gp = 0) {
    auto &grains = slots;
    auto &s = mixl;
    auto &s2 = mixr;
    std::size_t tt = grains[0].fs / dens;
    std::size_t vs = am.size();
    s.resize(vs);
    s2.resize(vs);
    if (st >= tt) {
        st -= tt;
        grains[num].trigger(gd, gp);
        num = num == slots.size() - 1 ? 0 : num + 1;
    }
    std::fill(s.begin(),s.end(),0);
    std::fill(s2.begin(),s2.end(),0);
    bool ch = 0;
    pan = (1. - pan)*.5f;
    S ppan = 1 - pan;
    std::size_t j;
    for (auto &grain: grains) {
        j = 0;
	grain.vsize(vs);
        for (auto &o : grain(am,f,pm)) {
          s[j] += o*ppan;
          s2[j++] += o*(1.-ppan);
	}
      ppan = ch ? pan : 1. - pan;
      ch = !ch ;
    }
    st+=vs;
    return s;
  }

  auto &channel(bool ch) {
    return ch ? mixr : mixl; 
  }

  
};
} // namespace Aurora
