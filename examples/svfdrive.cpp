// svfdrive.cpp:
// Nonlinear low pass sv filter example
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
#include "TwoPole.h"
#include <cstdlib>
#include <iostream>

double inline nlm(double s, double dr) { return std::tanh(s * dr) / dr; }

int main(int argc, const char *argv[]) {
  if (argc > 7) {
    double sr = argc > 8 ? std::atof(argv[8]) : Aurora::def_sr;
    auto dur = std::atof(argv[1]);
    auto a = std::atof(argv[2]);
    auto f = std::atof(argv[3]);
    auto cf = std::atof(argv[4]);
    auto rs = std::atof(argv[5]);
    auto drv = std::atof(argv[6]);
    auto typ = std::atof(argv[7]);

    Aurora::TableSet<double> wave(Aurora::SAW);
    Aurora::BlOsc<double> osc(&wave, sr);
    Aurora::TwoPole<double, nlm> fil(sr);
    double att = 0.1 * dur, dec = 0.2 * dur, sus = 0.7, rt = 0.1;
    auto func = Aurora::ads_gen<double>(att, dec, sus);
    Aurora::Env<double> env(func, rt, sr);
    bool gate = 1;
    auto d = 2 * (1 - (rs > 1 ? 1 : (rs < 0 ? 0 : rs)));
    for (int n = 0; n < osc.fs() * dur; n += osc.vsize())
      for (auto s : fil(osc(a, f), env(cf, 1000, gate), d, drv, typ)) {
        if (n > sr * (dur - rt))
          gate = 0;
        std::cout << s << std::endl;
      }
  } else
    std::cout << "usage: " << argv[0]
              << " dur(s) amp freq(Hz) cutoff(Hz) res drv typ [sr]"
              << std::endl;
  return 0;
}
