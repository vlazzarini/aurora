// drive.cpp:
// Nonlinear distortion example
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
#include "Func.h"
#include "Osc.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace Aurora;
inline float scl(float a, float b) { return a * b; }
struct Synth {
  float att, dec, sus;
  std::vector<float> wave;
  Env<float> env;
  Osc<float, lookupi> osc;
  Func<float, std::tanhf> drive;
  BinOp<float, scl> amp;

  Synth(float rt, float sr)
      : att(0.1f), dec(0.3f), sus(0.7f), wave(def_vsize),
        env(ads_gen(att, dec, sus), rt, sr), osc(&wave, sr), drive(), amp() {
    std::size_t n = 0;
    for (auto &s : wave) {
      s = sin<float>((1. / wave.size()) * n++);
    }
  };

  const std::vector<float> &operator()(float a, float f, float dr, bool gate) {
    return amp(env(0, a, gate), drive(osc(dr, f)));
  }
};

int main(int argc, const char *argv[]) {
  if (argc > 4) {
    double sr = argc > 5 ? std::atof(argv[5]) : def_sr;
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2]);
    auto f = std::atof(argv[3]);
    auto dr = std::atof(argv[3]);
    float rel = 0.1;
    Synth synth(rel, sr);
    bool gate = 1;
    for (int n = 0; n < sr * dur; n += def_vsize) {
      if (n > sr * (dur - rel))
        gate = 0;
      for (auto s : synth(a, f, dr, gate))
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) drive [sr]"
              << std::endl;
  return 0;
}
