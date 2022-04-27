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

template <typename S>
inline S sum(S a, S b) { return a + b; } 

template <typename S>
struct ConvReverb {
  IR<S>   ir1;
  IR<S>   ir2;  
  IR<S>   ir3;
  Conv<S> c1;
  Conv<S> c2;
  Conv<S> c3;
  BinOp<S, sum> mix;
 

  ConvReverb(const std::vector<S> &s1,
	     const std::vector<S> &s2,
	     const std::vector<S> &s3) :
    ir1(s1,32), ir2(s2,256), ir3(s3,4096),
    c1(&ir1), c2(&ir2), c3(&ir3) { }

  void reset(const std::vector<S> &s1,
	     const std::vector<S> &s2,
	     const std::vector<S> &s3) {
    ir1.reset(s1,32);
    c1.reset(&ir1);
    ir2.reset(s2,256);
    c2.reset(&ir2);
    ir3.reset(s3,4096);
    c3.reset(&ir3);
  }

  const std::vector<S> &operator()(const std::vector<S> &in, S g) {
    return mix(mix(c1(in,g*0.3),c2(in,g*0.3)),c3(in,g*0.3));
  }

};

template <typename S>
ConvReverb<S> create_reverb(std::vector<S> &imp) {
  if(imp.size() < 8192)
    imp.resize(8192);
  std::vector<S> s1(imp.begin()+32, imp.begin()+256);
  std::vector<S> s2(imp.begin()+256, imp.begin()+4096);
  std::vector<S> s3(imp.begin()+4096, imp.end());
  return ConvReverb(s1,s2,s3);
}

int main(int argc, const char **argv) {
  SF_INFO sfinfo, sfinfoir;
  SNDFILE *fpir, *fpin, *fpout;
  std::vector<double> impulse;
  sf_count_t n;

  if (argc > 4) {
    if ((fpir = sf_open(argv[1], SFM_READ, &sfinfoir)) != NULL) {
      if (sfinfoir.channels < 2) {
        impulse.resize(sfinfoir.frames);
        n = sf_read_double(fpir, impulse.data(), sfinfoir.frames);
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
        double g = atof(argv[4]);
        std::vector<double> buffer(def_vsize);
        ConvReverb<double> delay = create_reverb(impulse);
        Mix<double> mix;
        do {
          n = sf_read_double(fpin, buffer.data(), def_vsize);
          if (n) {
            buffer.resize(n);
            auto &out = mix(delay(buffer, g),buffer);
            sf_write_double(fpout, out.data(), n);
          } else
            break;
        } while (1);
        buffer.resize(def_vsize);
        n = impulse.size();
        std::fill(buffer.begin(), buffer.end(), 0.f);
        do {
          auto &out = delay(buffer, g);
          sf_write_double(fpout, out.data(), def_vsize);
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
