// follow.cpp:
// envelope follower example
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

#include "Func.h"
#include "OnePole.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;

template <typename S> S scl(S a, S b) { return a * b; }
template <typename S> S pos(S a) { return a < 0 ? -a : a; }

template <typename S> struct Follow {
  OnePole<S> filter;
  Func<S, pos<S>> abs;
  BinOp<S, scl<S>> amp;
  S fr;
  S g;

  Follow(S ga, S f = 10, S fs = def_sr, std::size_t vsize = def_vsize)
      : filter(fs, vsize), abs(vsize), amp(vsize), fr(f), g(ga){};

  const std::vector<S> &operator()(const std::vector<S> &in1,
                                   const std::vector<S> &in2) {
    return amp(g, amp(in1, filter(abs(in2), fr)));
  }
};

int main(int argc, const char **argv) {
  SF_INFO sfinfo1, sfinfo2;
  SNDFILE *fpin1, *fpin2, *fpout;
  sf_count_t n1, n2;

  if (argc > 4) {
    if ((fpin1 = sf_open(argv[1], SFM_READ, &sfinfo1)) == NULL) {
      std::cout << "could not open " << argv[1] << std::endl;
      return 1;
    }
    if ((fpin2 = sf_open(argv[2], SFM_READ, &sfinfo2)) == NULL) {
      std::cout << "could not open " << argv[1] << std::endl;
      sf_close(fpin1);
      return 1;
    }
    if (sfinfo1.samplerate != sfinfo2.samplerate) {
      std::cout << "sample rates do not match\n";
      sf_close(fpin1);
      sf_close(fpin2);
      return 1;
    }
    if (sfinfo2.channels > 1) {
      std::cout << "only mono files allowed\n";
      sf_close(fpin1);
      sf_close(fpin2);
      return 1;
    }

    fpout = sf_open(argv[3], SFM_WRITE, &sfinfo2);
    float gain = atof(argv[4]);
    std::vector<float> buffer1(def_vsize);
    std::vector<float> buffer2(def_vsize);
    Follow<float> follow(gain, 10.f, sfinfo2.samplerate);

    do {
      n1 = sf_read_float(fpin1, buffer1.data(), def_vsize);
      n2 = sf_read_float(fpin2, buffer2.data(), def_vsize);
      if (n1 && n2) {
        buffer1.resize(n1 < n2 ? n1 : n2);
        buffer2.resize(n1 < n2 ? n1 : n2);
        auto &out = follow(buffer1, buffer2);
        sf_write_float(fpout, out.data(), n1 < n2 ? n1 : n2);
      } else
        break;
    } while (1);
    sf_close(fpin1);
    sf_close(fpin2);
    sf_close(fpout);
    return 0;
  }
  std::cout << "usage: " << argv[0] << " infile envfile outfile gain \n"
            << std::endl;
  return -1;
}
