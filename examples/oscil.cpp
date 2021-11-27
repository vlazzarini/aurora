// oscil.cpp:
// Table lookup oscillator and envelope example
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
#include <cstdlib>
#include <iostream>

int main(int argc, const char *argv[]) {
  if (argc > 3) {
    double sr = argc > 4 ? std::atof(argv[4]) : Aurora::def_sr;
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2]);
    auto f = std::atof(argv[3]);
    std::vector<float> wave(16384);
    std::size_t n = 0;
    for (auto &s : wave)
      s = Aurora::sin<float>((1. / wave.size()) * n++);
    float att = 0.1f, dec = 0.3f, sus = 0.7f, rel = 0.1f;
    Aurora::Env<float> env(Aurora::ads_gen(att, dec, sus), rel, sr);
    Aurora::Osc<float> osc(Aurora::lookupi_gen(wave), sr);
    bool gate = 1;
    for (int n = 0; n < osc.fs() * dur; n += osc.vsize()) {
      if (n > osc.fs() * (dur - rel))
        gate = 0;
      for (auto s : env(osc(a, f), gate))
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) [sr]"
              << std::endl;
  return 0;
}
