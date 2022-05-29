// freqshift.cpp:
// Frequency shifter
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

#include "Quad.h"
#include "Osc.h"
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;
class FreqShifter {
  Quad<float> quad;
  Osc<float,phase> ph;
  std::vector<float> cost;
  std::vector<float> sint;  
  std::vector<float> up;
  std::vector<float> down;

public:
  FreqShifter(float sr, std::size_t vsize = def_vsize)
    :   quad(sr,vsize), ph(sr,vsize), cost(def_ftlen+1), sint(def_ftlen+1),
    up(vsize), down(vsize) {
    for(std::size_t n = 0; n < def_ftlen+1; n++) {
      cost[n] = std::cos(n*twopi/def_ftlen);
      sint[n] = std::sin(n*twopi/def_ftlen);
    }
  };


  void reset(float sr) {
    ph.reset(sr);
    quad.reset(sr);
  }


  const std::vector<float> &operator()(const std::vector<float> &in, float fr) {
    std::size_t n = 0;
    float re, im;
    up.resize(in.size());
    down.resize(in.size());
    ph.vsize(in.size());
    auto &phase = ph(1, fr);
    auto &inr = quad(in);
    auto &ini = quad.imag();
    for(auto &s : up) {
      re = inr[n]*lookupi(phase[n],&cost);
      im = ini[n]*lookupi(phase[n],&sint);
      s = re - im;
      down[n++] = re + im;
    }
    return up;
  }

  const std::vector<float> &downshift() {
    return down;
  }
  
};

int main(int argc, const char **argv) {
  SF_INFO sfinfo;
  SNDFILE *fpin, *fpout;
  if (argc > 3) {
    if ((fpin = sf_open(argv[1], SFM_READ, &sfinfo)) != NULL) {
      if (sfinfo.channels < 2) {
        std::size_t n = 0;
	float fr = std::atof(argv[3]);
        fpout = sf_open(argv[2], SFM_WRITE, &sfinfo);
        std::vector<float> in(def_vsize);
        FreqShifter shifter(sfinfo.samplerate);
        do {
          std::fill(in.begin(), in.end(), 0);
          n = sf_read_float(fpin, in.data(), def_vsize);
          auto &s = shifter(in, fr);
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

  std::cout << "usage: " << argv[0] << " infile outfile shift \n" << std::endl;
  return -1;
}
