// towave.c
// ASCII to RIFF-Wave conversion utility
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

#include <sndfile.h>
#include <stdio.h>
#include <stdlib.h>
#define SR 44100

int main(int argc, const char *argv[]) {
  SNDFILE *fp;
  SF_INFO sfinfo;
  const char *fname = "out.wav";
  float f[BUFSIZ];
  size_t n, i;

  sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
  sfinfo.samplerate = SR;
  sfinfo.channels = 1;

  switch (argc) {
  case 4:
    sfinfo.channels = atof(argv[3]);
  case 3:
    sfinfo.samplerate = atof(argv[2]);
  case 2:
    fname = argv[1];
  }

  fp = sf_open(fname, SFM_WRITE, &sfinfo);
  do {
    for (n = 0, i = 0; n < BUFSIZ; n++) {
      if (fscanf(stdin, "%f", &f[n]) > 0)
        ;
      else
        break;
    }
    sf_write_float(fp, f, n);
  } while (n);

  sf_close(fp);
  return 0;
}
