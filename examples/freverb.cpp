// freverb.cpp:
// reverb processing example
// depends on libsndfile
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
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;

template <typename S> inline S scl(S a, S b) { return a * b; }
template <typename S> struct Reverb {
static  constexpr S dt[4] = {0.037, 0.031, 0.029, 0.023};
static constexpr S adt[2] = {0.01, 0.0017};

  std::array<Del<S, lp_delay>, 4> combs;
  std::array<Del<S>, 2> apfs;
  Mix<S> mix;
  BinOp<S, scl> gain;
  std::array<std::vector<S>, 4> mem;
  std::array<S, 4> g;
  S rvt;

  void reverb_time(S rvt) {
    std::size_t n = 0;
    for (auto &gs : g)
      gs = std::pow(.001, dt[n++] / rvt);
  }

  void lp_freq(S lpf, S fs) {
    for (auto &m : mem) {
      m[0] = 0;
      double c = 2. - std::cos(2 * M_PI * lpf / fs);
      m[1] = sqrt(c * c - 1.f) - c;
    }
  }

  void reset(S rvt, S lpf, S fs) {
    std::size_t n = 0;
    for (auto &obj : combs)
      obj.reset(dt[n++], fs);
    apfs[0].reset(adt[0], fs);
    apfs[1].reset(adt[1], fs);
    reverb_time(rvt);
    lp_freq(lpf, fs);
  }

  Reverb(S rvt, S lpf, S fs = def_sr, std::size_t vsize = def_vsize)
      : combs({Del<S, lp_delay>(dt[0], fs, vsize),
               Del<S, lp_delay>(dt[1], fs, vsize),
               Del<S, lp_delay>(dt[2], fs, vsize),
               Del<S, lp_delay>(dt[3], fs, vsize)}),
        apfs({Del<S>(adt[0], fs, vsize), Del<S>(adt[1], fs, vsize)}),
        mix(vsize), gain(vsize), mem({std::vector<S>(2), std::vector<S>(2),
                                      std::vector<S>(2), std::vector<S>(2)}),
        g({0, 0, 0, 0}) {
    reverb_time(rvt);
    lp_freq(lpf, fs);
  };

  const std::vector<S> &operator()(const std::vector<S> &in, S rmx) {
    S ga0 = 0.7;
    S ga1 = 0.7;
    auto &s = gain(0.25, mix(combs[0](in, 0, g[0], 0, &mem[0]),
                             combs[1](in, 0, g[1], 0, &mem[1]),
                             combs[2](in, 0, g[2], 0, &mem[2]),
                             combs[3](in, 0, g[3], 0, &mem[3])));
    return mix(in, gain(rmx, apfs[1](apfs[0](s, 0, ga0, -ga0), 0, ga1, -ga1)));
  }
};

int main(int argc, const char **argv) {
  SF_INFO sfinfo;
  SNDFILE *fpin, *fpout;
  int n;

  if (argc > 5) {
    if ((fpin = sf_open(argv[1], SFM_READ, &sfinfo)) != NULL) {
      if (sfinfo.channels < 2) {
        fpout = sf_open(argv[2], SFM_WRITE, &sfinfo);
        float rvt = atof(argv[3]);
        float rmx = atof(argv[4]);
        float cf = atof(argv[5]);
        std::vector<float> buffer(def_vsize);
        Reverb<float> reverb(rvt, cf, sfinfo.samplerate);
        do {
          std::fill(buffer.begin(), buffer.end(), 0);
          n = sf_read_float(fpin, buffer.data(), def_vsize);
          if (n) {
            buffer.resize(n);
            auto &out = reverb(buffer, rmx);
            sf_write_float(fpout, out.data(), n);
          } else
            break;
        } while (1);
        std::cout << buffer.size() << std::endl;
        buffer.resize(def_vsize);
        std::cout << buffer.size() << std::endl;
        n = sfinfo.samplerate * rvt;
        std::fill(buffer.begin(), buffer.end(), 0);
        do {
          auto &out = reverb(buffer, rmx);
          sf_write_float(fpout, out.data(), def_vsize);
          n -= def_vsize;
        } while (n > 0);
        sf_close(fpin);
        sf_close(fpout);
        return 0;
      } else
        std::cout << "only mono soundfiles permitted\n";
      sf_close(fpin);
      return 1;
    } else
      std::cout << "could not open " << argv[1] << std::endl;
    return 1;
  }
  std::cout << "usage: " << argv[0]
            << " infile outfile reverb_time reverb_amount lpf \n"
            << std::endl;
  return -1;
}
