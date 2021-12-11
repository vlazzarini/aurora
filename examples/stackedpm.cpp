// stackedpm.cpp:
// Stacked PM example
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

namespace Aurora {

inline double scl(double a, double b) { return a * b; }
template <typename S> struct StackedPM {
  Osc<S, sin<S>> mod0;
  Osc<S, sin<S>> mod1;
  Osc<S> car;
  S att, dec, sus;
  Env<S> env;
  BinOp<S, scl> amp;
  S o2pi;

  StackedPM(S rel, S fs = (S)def_sr, std::size_t vsize = def_vsize)
      : mod0(fs, vsize), mod1(fs, vsize), car(fs, vsize), att(0), dec(0),
        sus(0), env(ads_gen(att, dec, sus), rel, fs), amp(), o2pi(1. / twopi){};

  void release(S rel) { env.release(rel); }

  std::size_t vsize() { return car.vsize(); }

  S fs() { return car.fs(); }

  const std::vector<S> &operator()(S a, S fc, S fm0, S fm1, S z0, S z1,
                                   bool gate, std::size_t vsiz = 0) {
    if (vsiz) {
      env.vsize(vsiz);
      mod0.vsize(vsiz);
    }
    auto &e = env(gate);
    auto ke = e[0];
    auto &s0 = mod0(z0 * o2pi * (ke + 1), fm0);
    auto &s1 = mod1(z1 * o2pi * ke, fm1, s0);
    return amp(e, car(a, fc, s1));
  }
};
} // namespace Aurora

int main(int argc, const char *argv[]) {
  if (argc > 3) {
    double sr = argc > 4 ? std::atof(argv[4]) : Aurora::def_sr;
    auto dur = std::atof(argv[1]);
    auto amp = std::atof(argv[2]);
    auto fr = std::atof(argv[3]);
    Aurora::StackedPM<double> fm(0.1, sr);
    fm.att = 0.05;
    fm.dec = 0.8;
    fm.sus = 0.6;
    bool gate = 1;
    for (int n = 0; n < fm.fs() * (dur + 0.91); n += fm.vsize()) {
      if (n > dur * fm.fs())
        gate = 0;
      const std::vector<double> &out =
          fm(amp, fr, fr * .999f, fr * 1.001f, 2, 3, gate);
      for (auto s : out)
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) [sr]"
              << std::endl;
  return 0;
}
