// reverb.cpp:
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

#include "Conv.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;

int main(int argc, const char **argv) {
  SF_INFO sfinfo, sfinfoir;
  SNDFILE *fpir, *fpin, *fpout;
  std::vector<float> impulse;
  sf_count_t n;

  if (argc > 4) {
    if ((fpir = sf_open(argv[1], SFM_READ, &sfinfoir)) != NULL) {
      if (sfinfoir.channels < 2) {
        impulse.resize(sfinfoir.frames);
        n = sf_read_float(fpir, impulse.data(), sfinfoir.frames);
        if (!n) {
          std::cout << "error reading " << argv[1] << std::endl;
          sf_close(fpir);
          return 1;
        }
      } else {
        std::cout << "only mono soundfiles permitted\n";
        return 1;
      }
    } else {
      std::cout << "error opening" << argv[1] << std::endl;
      return 1;
    }

    sf_close(fpir);
    if ((fpin = sf_open(argv[2], SFM_READ, &sfinfo)) != NULL) {

      if (sfinfoir.samplerate != sfinfo.samplerate) {
        std::cout << "sample rates do not match\n";
        sf_close(fpin);
        return 1;
      }

      if (sfinfo.channels < 2) {
        fpout = sf_open(argv[3], SFM_WRITE, &sfinfo);
        float g = atof(argv[4]);
        std::vector<float> buffer(def_vsize);
        IR<float> ir(impulse, 512);
        Conv<float> delay(&ir);
        Mix<float> mix;
        do {
          n = sf_read_float(fpin, buffer.data(), def_vsize);
          if (n) {
            buffer.resize(n);
            auto &out = mix(delay(buffer, g), buffer);
            sf_write_float(fpout, out.data(), n);
          } else
            break;
        } while (1);
        buffer.resize(def_vsize);
        n = impulse.size();
        std::fill(buffer.begin(), buffer.end(), 0.f);
        do {
          auto &out = delay(buffer, g);
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
      std::cout << "could not open " << argv[2] << std::endl;
    return 1;
  }
  std::cout << "usage: " << argv[0] << " irfile infile outfile rev_gain \n"
            << std::endl;
  return -1;
}
