
//
// (c) V Lazzarini, 2025
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
#include "Func.h"
#include "BlOsc.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace Aurora;

struct Synth {
  inline static const float smax = 8.f;
  inline static std::vector<float> sigmoid{std::vector<float>(0)};
  inline static TableSet<float> wave{SAW};


  static float sat(float a) {
    return cubic_interp_lim((a / smax + .5) * sigmoid.size(), sigmoid);
  }
  static float scl(float a, float b) { return a * b; }

  float att, dec, sus;

  Env<float> env;
  BlOsc<float, lookupi> osc;
  Func<float, sat> drive;
  BinOp<float, scl> amp;

  Synth(float rt, float sr)
      : att(0.1f), dec(0.3f), sus(0.7f),
        env(ads_gen(att, dec, sus), rt, sr / def_vsize, 1), osc(&wave, sr),
        drive(), amp() {
    if (sigmoid.size() == 0) {
      std::size_t n = 0;
      sigmoid.resize(def_ftlen);
      for (auto &s : sigmoid) {
        s = std::tanhf((smax / sigmoid.size()) * n++ - smax / 2);
      }
    }
  }

  const std::vector<float> &operator()(float a, float f, float dr, bool gate,
                                       std::size_t vsiz = 0) {
    if (vsiz) {
      osc.vsize(vsiz);
      env.vsize(vsiz);
    }
    return amp(a * env(gate)[0], drive(osc(dr, f)));
  }
};

int main(int argc, const char *argv[]) {
  if (argc > 4) {
    double sr = argc > 6 ? std::atof(argv[6]) : def_sr;
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2])*0.5;
    auto f1 = std::atof(argv[3]);
    auto f2 = std::atof(argv[4]);
    auto dr = std::atof(argv[5]);
    float rel = 0.1;

    std::vector<Synth> sy;
    std::vector<float> sig(def_vsize);
    sy.push_back(Synth(rel, sr));
    sy.push_back(Synth(rel, sr));
    bool gate = 1;
    for (int n = 0; n < sr * dur; n += def_vsize) {
      int i = 0;
      std::fill(sig.begin(), sig.end(), 0.f);
      if (n > sr * (dur - rel))
        gate = 0;
      auto &s1 = sy[0](a, f1, dr, gate);
      for (auto s : sy[1](a, f2, dr, gate)) 
        std::cout << s + s1[i++] << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq1(Hz) freq2(Hz) drive [sr]"
              << std::endl;
  return 0;
}
