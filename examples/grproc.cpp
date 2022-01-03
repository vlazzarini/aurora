// grproc.cpp:
// grain processing example
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

#include "grain.h"
#include <cstdlib>
#include <iostream>
#include <sndfile.h>
#include <vector>

using namespace Aurora;

int main(int argc, const char **argv) {
  SF_INFO sfinfo_in, sfinfo;
  SNDFILE *fpin, *fpout;

  if (argc > 7) {
    if ((fpin = sf_open(argv[1], SFM_READ, &sfinfo_in)) != NULL) {
      if (sfinfo_in.channels < 2) {
        sfinfo = sfinfo_in;
        fpout = sf_open(argv[2], SFM_WRITE, &sfinfo);
        float a = atof(argv[3]);
        float p = atof(argv[4]);
        float t = atof(argv[5]) * def_vsize / double(sfinfo_in.samplerate);
        float gdur = atof(argv[6]);
        float ol = atof(argv[7]) ;
        int dm = argc > 8 ? atoi(argv[8]) : def_vsize;
        double ts = t >= 0 ? 0.f : sfinfo_in.frames / double(sfinfo_in.samplerate);
        std::vector<float> wave(sfinfo_in.frames);
        sf_read_float(fpin, wave.data(), sfinfo_in.frames);
        sf_close(fpin);
        GrainGen<float> grain(wave, std::ceil(ol), sfinfo_in.samplerate, dm);
        do {
          auto &out = grain(a, p, ol / gdur, gdur, ts);
          ts += t; 
          sf_write_float(fpout, out.data(), def_vsize);
        } while (t < 0 ? ts >= 0
                 : ts < sfinfo_in.frames / double(sfinfo_in.samplerate));
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
            << " infile outfile amp pitchscale timescale grainsize(s) overlap "
               "[decimation]\n"
            << std::endl;
  return -1;
}
