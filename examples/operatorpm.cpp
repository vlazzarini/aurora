// operatorpm.cpp:
// Operator PM example
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

using namespace Aurora;
template <typename S> struct Opm {
  static S scl(S a, S b) { return a * b; }
  Osc<S, sin<S>> osc;
  Env<S> env;
  BinOp<S, scl> amp;
  S att, dec, sus;
  S o2pi;

  Opm(S rel, S fs = def_sr, std::size_t vsize = def_vsize)
      : osc(fs), env(att, dec, sus, rel, fs), att(0), dec(0), sus(0),
        o2pi(1 / twopi) {}

  void release(S rel) { env.release(rel); }

  auto &operator()(S a, S f, bool gate, std::size_t vsiz = 0) {
    if (vsiz)
      osc.vsize(vsiz);
    return env(osc(a, f), gate);
  }

  auto &operator()(S a, S f, const std::vector<S> &pm, bool gate) {
    return env(osc(a, f, amp(o2pi, pm)), gate);
  }
};

int main(int argc, const char *argv[]) {
  if (argc > 3) {
    double sr = argc > 4 ? std::atof(argv[4]) : Aurora::def_sr;
    auto dur = std::atof(argv[1]);
    auto amp = std::atof(argv[2]);
    auto fr = std::atof(argv[3]);
    Opm<double> op1(0.1, sr);
    op1.att = 0.5;
    op1.dec = 0.8;
    op1.sus = 0.3;
    Opm<double> op2(0.1, sr);
    op2.att = 0.001;
    op2.dec = 0.1;
    op2.sus = 0.7;
    bool gate = 1;
    for (int n = 0; n < sr * (dur + 0.91); n += def_vsize) {
      if (n > dur * sr)
        gate = 0;
      auto &out = op2(amp, fr, op1(6, 2 * fr, gate), gate);
      for (auto s : out)
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) [sr]"
              << std::endl;
  return 0;
}
