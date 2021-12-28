// noise.cpp:
// White noise and envelope example
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
#include <cstdlib>
#include <iostream>

using namespace Aurora;
struct Synth {
  static float randf(float a) {
    return a * (((std::rand()) / float(RAND_MAX) * 2) - 1.f);
  }

  float att, dec, sus;
  Env<float> env;
  Func<float, randf> noise;

  Synth(float rt, float sr)
      : att(0.f), dec(0.f), sus(0.f), env(ads_gen(att, dec, sus), rt, sr),
        noise(){};

  const std::vector<float> &operator()(float a, bool gate,
                                       std::size_t vsiz = 0) {
    if (vsiz)
      noise.vsize(vsiz);
    return env(noise(a), gate);
  }
};

int main(int argc, const char *argv[]) {
  if (argc > 2) {
    double sr = argc > 3 ? std::atof(argv[3]) : def_sr;
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2]);
    float rel = 0.1;
    Synth synth(rel, sr);
    synth.att = .1f;
    synth.dec = .3f;
    synth.sus = .7f;
    bool gate = 1;
    for (int n = 0; n < sr * (dur + rel); n += def_vsize) {
      if (n > sr * dur)
        gate = 0;
      for (auto s : synth(a, gate))
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp [sr]" << std::endl;
  return 0;
}
