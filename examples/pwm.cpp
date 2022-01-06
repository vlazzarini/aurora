// pwm.cpp:
// variable pulse width example
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

#include "FourPole.h"
#include "BlOsc.h"
#include "Env.h"
#include <cstdlib>
#include <iostream>

using namespace Aurora;
float add(float a, float b) { return a + b; }
struct Synth {
  float att, dec, sus;
  TableSet<float> wave;
  Env<float> env;
  BlOsc<float, lookupi<float>> osc1;
  BlOsc<float, lookupi<float>> osc2;
  BinOp<float, add> mix;

  Synth(float rt, float sr)
      : att(0.f), dec(0.f), sus(0.f), wave(SAW),
        env(ads_gen(att, dec, sus), rt, sr), osc1(&wave, sr), osc2(&wave, sr),
        mix(){};

  const std::vector<float> &operator()(float a, float f, float pwm, bool gate,
                                       std::size_t vsiz = 0) {
    if (vsiz) {
      osc1.vsize(vsiz);
      osc2.vsize(vsiz);
    }
    float off = a * (2 * pwm - 1.f);
    auto &m = mix(mix(osc1(a, f, pwm), osc2(-a, f)), off);
    return env(m, gate);
  }
};

int main(int argc, const char *argv[]) {
  if (argc > 4) {
    double sr = argc > 5 ? std::atof(argv[5]) : def_sr;
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2]);
    auto f = std::atof(argv[3]);

    auto pwm = std::atof(argv[4]);
    float rel = 0.1;
    Synth synth(rel, sr);
    Aurora::FourPole<float> fil(sr);
    float att = 0.01 * dur, dec = 0.1 * dur, sus = 0.1, rt = 0.1;
    Aurora::Env<float> env(att, dec, sus, rt, sr);
    Aurora::Osc<float> lfo(sr);
    synth.att = .1f;
    synth.dec = .3f;
    synth.sus = .7f;
    bool gate = 1;
    float res = 0.5;
    float cf = 10000.f;
    for (int n = 0; n < sr * (dur + rel); n += def_vsize) {
      if (n > sr * dur)
        gate = 0;
      for (auto s : fil(synth(a, f, lfo(pwm/2,0.9)[0] + 0.5 , gate), env(f, cf - f, gate), res))
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) pwm(0-1) [sr]"
              << std::endl;
  return 0;
}
