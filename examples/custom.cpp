// custom.cpp:
// Custom bandlimited waveform generation example
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

#include "BlOsc.h"
#include "Env.h"
#include "flute.h"
#include <cstdlib>
#include <iostream>

int main(int argc, const char *argv[]) {
  if (argc > 3) {
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2]);
    auto f = std::atof(argv[3]) * flute::ratio;
    Aurora::TableSet<float> wave(flute::wave, flute::base, flute::fs);
    Aurora::BlOsc<float> osc(&wave, flute::fs);
    std::vector<float> brkpts({0.1f, 1.f, 0.5f, 0.1f, 1.f, 0.7f});
    Aurora::Env<float> env(Aurora::env_gen<float>(brkpts), 0.1f,
                           flute::fs / Aurora::def_vsize, 1);
    bool gate = 1;
    for (int n = 0; n < osc.fs() * (dur + 0.1f); n += osc.vsize()) {
      if (n > osc.fs() * dur)
        gate = 0;
      for (auto s : osc(env(0, a, gate)[0], f))
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz)" << std::endl;
  return 0;
}
