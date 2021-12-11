// stackedfm.cpp:
// Stacked FM example
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

#include "Osc.h"
#include <cstdlib>
#include <iostream>

namespace Aurora {

inline double scl(double a, double b) { return a * b; }
inline double mix(double a, double b) { return a + b; }

template <typename S> class StackedFM {
  Osc<S> mod0;
  Osc<S> mod1;
  Osc<S> car;
  BinOp<S, scl> amp;
  BinOp<S, mix> add;

public:
  StackedFM(S fs = (S)def_sr, std::size_t vsize = def_vsize)
      : mod0(fs, vsize), mod1(fs, vsize), car(fs, vsize), amp(vsize),
        add(vsize){};

  std::size_t vsize() { return car.vsize(); }

  S fs() { return car.fs(); }

  const std::vector<S> &operator()(S a, S fc, S fm0, S fm1, S z0, S z1,
                                   std::size_t vsiz = 0) {
    if (vsiz)
      mod0.vsize(vsiz);
    auto &s0 = add(fm1, mod0(z0 * fm0, fm0));
    auto &s1 = add(fc, mod1(amp(z1, s0), s0));
    return car(a, s1);
  }
};
} // namespace Aurora

int main(int argc, const char *argv[]) {
  if (argc > 3) {
    double sr = argc > 4 ? std::atof(argv[4]) : Aurora::def_sr;
    auto dur = std::atof(argv[1]);
    auto amp = std::atof(argv[2]);
    auto fr = std::atof(argv[3]);
    Aurora::StackedFM<double> fm(sr);
    for (int n = 0; n < fm.fs() * dur; n += fm.vsize()) {
      const std::vector<double> &out = fm(amp, fr, fr, fr, 3, 2);
      for (auto s : out)
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) [sr]"
              << std::endl;
  return 0;
}
