// equaliser.cpp:
// equaliser example
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

#include "Eq.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;
template <typename S> struct Equaliser {
  std::vector<Eq<S>> eq;
  Equaliser(std::size_t nf, S fs = def_sr, std::size_t vsize = def_vsize)
      : eq(nf, Eq<S>(fs, vsize)){};

  const std::vector<S> &operator()(const std::vector<S> &in,
                                   const std::vector<S> &gs,
                                   const std::vector<S> &f,
                                   const std::vector<S> &bw) {
    auto s = &in;
    for (std::size_t n = 0; n < gs.size(); n++)
      s = &eq[n](*s, gs[n], f[n], bw[n]);
    return *s;
  }
};

int main(int argc, const char **argv) {
  SF_INFO sfinfo;
  SNDFILE *fpin, *fpout;
  FILE *fp;
  sf_count_t n;

  if (argc > 3) {
    if ((fp = fopen(argv[1], "r")) == NULL) {
      std::cout << "could not open " << argv[1] << std::endl;
      return 1;
    }
    if ((fpin = sf_open(argv[2], SFM_READ, &sfinfo)) == NULL) {
      std::cout << "could not open " << argv[1] << std::endl;
      fclose(fp);
      return 1;
    }

    if (sfinfo.channels > 1) {
      std::cout << "only mono files allowed\n";
      fclose(fp);
      sf_close(fpin);
      return 1;
    }

    fpout = sf_open(argv[3], SFM_WRITE, &sfinfo);
    std::vector<float> f;
    std::vector<float> b;
    std::vector<float> g;
    int i = 0;
    do {
      float tg, tf, tb; 
      n = fscanf(fp, "%f  %f  %f\n", &tg, &tf, &tb);
      if (n > 0) {
        g.push_back(tg);
        f.push_back(tf);
        b.push_back(tb);
        printf("band %d - g:%.3f  cf:%.1f  bw:%.1f\n", ++i, tg, tf, tb);
      }
    } while (n > 0);
    fclose(fp);

    std::vector<float> buffer(def_vsize);
    Equaliser<float> eq(g.size(), sfinfo.samplerate);

    do {
      n = sf_read_float(fpin, buffer.data(), def_vsize);
      if (n) {
        buffer.resize(n);
        auto &out = eq(buffer, g, f, b);
        sf_write_float(fpout, out.data(), n);
      } else
        break;
    } while (1);
    sf_close(fpin);
    sf_close(fpout);
    return 0;
  }
  std::cout << "usage: " << argv[0] << " paramfile infile outfile \n"
            << std::endl;
  return -1;
}
