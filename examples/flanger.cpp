// filter.cpp:
// Low pass resonant filter processing example
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
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;

int main(int argc, const char **argv) {
  SF_INFO sfinfo;
  SNDFILE *fpin, *fpout;
  int n;

  if (argc > 6) {
    if ((fpin = sf_open(argv[1], SFM_READ, &sfinfo)) != NULL) {
      if (sfinfo.channels < 2) {
        fpout = sf_open(argv[2], SFM_WRITE, &sfinfo);
        float mxdel = atof(argv[3]) / 1000.;
        mxdel = mxdel < 10.f ? mxdel : 10.f;
        float fr = atof(argv[4]);
        float fdb = atof(argv[5]);
        float g = atof(argv[6]);
        std::vector<float> buffer(def_vsize);
        auto lfofun = [](float x) -> float {
          return Aurora::cos<float>(x) * 0.49 + 0.51;
        };
        Osc<float> lfo(lfofun, sfinfo.samplerate);
        auto delf = vdelayc<float>;
        Del<float> delay(0.01, delf, sfinfo.samplerate);
        BinOp<float> gain([](float a, float b) -> float { return a * b; });
        do {
          std::fill(buffer.begin(), buffer.end(), 0);
          n = sf_read_float(fpin, buffer.data(), def_vsize);
          auto &dlt = lfo(mxdel, fr);
          auto &out = gain(delay(buffer, dlt, fdb), g);
          sf_write_float(fpout, out.data(), n);
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
  std::cout << "usage: " << argv[0]
            << " infile outfile maxdel(ms) lfofr(Hz) fdb gain \n"
            << std::endl;
  return -1;
}
