// Chorus.cpp:
// Chorus effect
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
#include "Osc.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;
struct StereoChorus {
  static float ofs(float a, float b) { return a + b; }
  std::vector<float> out;
  std::array<Osc<float>, 2> lfo;
  std::array<Del<float, vdelayi>, 2> delay;
  BinOp<float, ofs> offs;
  

  StereoChorus(float sr, std::size_t vsize = def_vsize)
    :    out(vsize*2),
        lfo{Osc<float>(sr, vsize), Osc<float>(sr, vsize)},
        delay{Del<float, vdelayi>(0.1, sr, vsize),
              Del<float, vdelayi>(0.1, sr, vsize)},
        offs(vsize){};

  auto &chorus(const std::vector<float> &in, float fr, float d, int chn) {
    lfo[chn].vsize(in.size());
    return delay[chn](in, offs(d, lfo[chn](d * 0.1, fr)), 0, 1);
  }

  auto &operator()(const std::vector<float> &in, float width) {
          bool chn = 0;
          std::size_t n = 0;
          out.resize(in.size()*2);  
          auto &l = chorus(in, .93f, .017f, 0);
          auto &r = chorus(in, .87f, .013f, 1);
          for (auto &b : out) {
            b = (chn = !chn) ? l[n] * width - r[n] * (1-width) 
                             : r[n] * width - l[n] * (1-width);
            n = !chn ? n + 1 : n;
          }
          return out;
  }
};

int main(int argc, const char **argv) {
  SF_INFO sfinfo;
  SNDFILE *fpin, *fpout;
  if (argc > 2) {
    if ((fpin = sf_open(argv[1], SFM_READ, &sfinfo)) != NULL) {
      if (sfinfo.channels < 2) {
        std::size_t n = 0;
        sfinfo.channels = 2;
        fpout = sf_open(argv[2], SFM_WRITE, &sfinfo);
        std::vector<float> in(def_vsize);
        StereoChorus chorus(sfinfo.samplerate);
        do {
          std::fill(in.begin(), in.end(), 0);
          n = sf_read_float(fpin, in.data(), def_vsize);
          auto &s = chorus(in, 0.7f);
          sf_writef_float(fpout, s.data(), n);
        } while (n);
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

  std::cout << "usage: " << argv[0] << " infile outfile \n" << std::endl;
  return -1;
}
