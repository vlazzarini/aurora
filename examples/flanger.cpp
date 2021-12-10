// flanger.cpp:
// Flanger effect
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

inline float scl(float a, float b) { return a * b; }
inline float lfofun(double x, const std::vector<float> *nop) {
  return cos<float>(x) * 0.46 + 0.54;
}

struct Flanger {
  Osc<float, lfofun> lfo;
  Del<float, vdelayi> delay;
  BinOp<float, scl> gain;
  float mxdel;

  Flanger(float maxdt, float sr)
      : lfo(sr), delay(maxdt, sr), gain(), mxdel(maxdt) {}

  const std::vector<float> &operator()(const std::vector<float> &in, float fr,
                                       float fdb, float g) {
    lfo.vsize(in.size());
    return gain(delay(in, lfo(mxdel, fr), fdb), g);
  }
};

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
        Flanger flanger(mxdel, sfinfo.samplerate);
        do {
          std::fill(buffer.begin(), buffer.end(), 0);
          n = sf_read_float(fpin, buffer.data(), def_vsize);
          auto &out = flanger(buffer, fr, fdb, g);
          sf_write_float(fpout, out.data(), n);

        } while (n);
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
            << " infile outfile maxdel(ms) lfofr(Hz) fdb gain \n"
            << std::endl;
  return -1;
}
