// karplus.cpp:
// Karplus-Strong synthesis example
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

#include "Del.h"
#include "Env.h"
#include <cstdlib>
#include <functional>
#include <iostream>

namespace Aurora {

template <typename S>
inline S kp(S rp, std::size_t wp, const std::vector<S> &del,
            std::vector<S> *mem) {
  std::size_t ds = del.size();
  S x = linear_interp(rpos(rp, wp, ds), del);
  if (mem != nullptr) {
    auto &dd = *mem;
    S y = (x + dd[0]) * 0.5;
    dd[0] = x;
    return y;
  } else
    return x;
}

template <typename S> inline S gat(S a, S d, S s, double t, S e, S ts) {
  return 1;
}

template <typename S> inline S scl(S a, S b) { return a * b; }

template <typename S> struct Karplus {

  Del<S, kp> delay;
  Env<S, gat> env;
  BinOp<S, scl> amp;
  std::vector<S> mem;
  std::vector<S> rnd;
  std::vector<S> in;
  S sr;
  bool gate;
  S g, ff, ddt;
  S twopiosr;

  void fill_delay() {
    delay(rnd, 0);
    mem[0] = 0;
  }

  Karplus(S fs = def_sr, std::size_t vsiz = def_vsize)
      : delay(0.05, fs, vsiz), env(nullptr, 0.1, fs, vsiz), amp(vsiz), mem(1),
        rnd(fs * 0.05), in(vsiz), sr(fs), gate(0), g(1), ff(0), ddt(0),
        twopiosr(2 * M_PI / sr) {
    for (auto &s : rnd) {
      s = 2 * (std::rand() / (S)RAND_MAX) - 1;
    }
  }

  void reset(S fs) {
    sr = fs;
    rnd.resize(fs * 0.05);
    for (auto &s : rnd) {
      s = 2 * (std::rand() / (S)RAND_MAX) - 1;
    }
    twopiosr = 2 * M_PI / sr;
    delay.reset(0.05, fs);
    ff = 0;
    ddt = 0;
    g = 1;
  }

  void release(S rel) { env.release(rel); }

  std::size_t vsize() { return in.size(); }

  void vsize(std::size_t n) { in.resize(n); }

  S fs() { return sr; }

  void decay(S fr, S dt) {
    S gf = std::pow(10, -60. / (20 * fr * dt));
    S gg = std::cos(fr * twopiosr);
    if (gf < gg)
      g = gf / gg;
    else
      g = 1.;
    ff = fr;
    ddt = dt;
  }

  void note_on() {
    fill_delay();
    gate = 1;
  }

  void note_off() { gate = 0; }

  const std::vector<S> &operator()(S a, S fr, S dt, std::size_t vsiz = 0) {
    if (vsiz)
      vsize(vsiz);
    fr = fr > 20 ? fr : 20;
    if (ff != fr || ddt != dt)
      decay(fr, dt);
    return amp(a, env(delay(in, 1 / fr - 1 / (2 * sr), g, 0, &mem), gate));
  }
};
} // namespace Aurora

int main(int argc, const char *argv[]) {
  if (argc > 3) {
    double sr = argc > 4 ? std::atof(argv[4]) : Aurora::def_sr;
    auto dec = std::atof(argv[1]);
    auto amp = std::atof(argv[2]);
    auto fr = std::atof(argv[3]);
    Aurora::Karplus<double> pluck(sr);
    pluck.note_on();
    for (int n = 0; n < pluck.fs() * (dec + 0.1); n += pluck.vsize()) {
      if (n > pluck.fs() * dec)
        pluck.note_off();
      auto &out = pluck(amp, fr, dec);
      for (auto s : out)
        std::cout << s << std::endl;
    }
  } else
    std::cout << "usage: " << argv[0] << " dur(s) amp freq(Hz) [sr]"
              << std::endl;
  return 0;
}
